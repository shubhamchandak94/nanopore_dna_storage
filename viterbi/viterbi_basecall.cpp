#include <omp.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <map>

const uint8_t NBASE = 4;
const char int2base[NBASE] = {'A', 'C', 'G', 'T'};
const bool base2bit[NBASE][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
const uint8_t nstate_crf = 8;
// states are A+,C+,G+,T+,A-,C-,G-,T- (+ = flip, - = flop), you can enter from a
// different base only into a flip base.

typedef std::array<std::array<float, nstate_crf>, NBASE + 1> crf_mat_t;
// for each input state, we have five (4+1) possible output states: e.g., A+ ->
// A+,C+,G+,T+,A-; A- -> A+,C+,G+,T+,A-; A-; etc.

struct path_t {
  std::string msg;
  uint32_t st_pos;
  uint32_t st_conv;
  uint8_t st_crf;
  float score;
};
// structure storing path for list decoding
// note that st_pos and st_conv are not really needed (can be obtained from msg)
// but we keep them for convenience
// msg stored as string to allow hashing

// for parallel LVA
const uint8_t BITSET_SIZE = 192;  // 24 bytes
typedef std::bitset<BITSET_SIZE> bitset_t;
struct LVA_path_t {
  bitset_t msg;
  float score;
  bool operator==(const LVA_path_t &other) const {
     return msg == other.msg;
  }
};

// convolutional code related parameters
const uint8_t mem_conv =
//     6; // mem 6 from CCSDS
//    8;  // mem 8 from GL paper
    11;  // mem 11 from GL paper
// 14; // mem 14 from GL paper
const uint32_t nstate_conv = 1 << mem_conv;
const uint8_t n_out_conv = 2;
typedef std::array<std::array<uint32_t, 2>, nstate_conv> conv_arr_t;
const uint32_t G[n_out_conv] =  // octal
//                                {0171, 0133};  // mem 6 from CCSDS
//    {0515, 0677};               // mem 8 from GL paper
    {05537, 06131};  // mem 11 from GL paper
// {075063, 056711}; // mem 14 from GL paper
const uint32_t initial_state_conv =  // binary
                                     //    0;                               // 0
                                     //    initial state
//                                          0b100101; // mem 6
//    0b10010110;                      // mem 8
    0b10010110001;  // mem 11
// 0b10010110001101; // mem 14

// when using sync_markers, initial_state = 0 should work just fine
const uint32_t sync_marker_length = 0;
// 1
// 2;
// 3;
const char sync_marker[sync_marker_length] = {};
// {1};
// {1,0};
// {1, 1, 0};
const uint32_t sync_marker_period = 9;

void generate_conv_arrays(conv_arr_t &prev_state, conv_arr_t &next_state,
                          conv_arr_t *output);

std::vector<bool> encode(const std::vector<bool> &msg,
                         const conv_arr_t &next_state,
                         const conv_arr_t *output);

std::vector<bool> read_bit_array(const std::string &infile);

void write_bit_array(const std::vector<bool> &outvec,
                     const std::string &outfile);

void write_char_array(const std::vector<char> &vec, const std::string &outfile);

void write_bit_array_in_bases(const std::vector<bool> &outvec,
                              const std::string &outfile);

template <class T>
void write_vector(const std::vector<T> &outvec, const std::string &outfile) {
  // write values in vector, one per line
  std::ofstream fout(outfile);
  for (auto v : outvec) {
    fout << v << "\n";
  }
  fout.close();
}

float logsumexpf(float x, float y) {
  float max_x_y = std::max(x, y);
  return max_x_y + logf(expf(x - max_x_y) + expf(y - max_x_y));
}

std::vector<crf_mat_t> read_crf_post(const std::string &infile);

std::vector<std::vector<uint8_t>> read_vocab_file(const std::string &infile);

uint8_t to_idx_crf_in_post(uint8_t st2_crf);

std::vector<char> decode_post_no_conv(const std::vector<crf_mat_t> &post);

std::vector<char> crfpath_to_basecall(const std::vector<uint8_t> &path);

uint32_t get_state_idx(const uint32_t st_pos, const uint32_t st_conv,
                       const uint32_t st_crf);

uint32_t state_idx_to_pos(const uint32_t st);

std::vector<bool> decode_post_conv(const std::vector<crf_mat_t> &post,
                                   const conv_arr_t &prev_state,
                                   const conv_arr_t &next_state,
                                   const conv_arr_t *output,
                                   const uint32_t msg_len);

std::vector<path_t> decode_post_conv_list(const std::vector<crf_mat_t> &post,
                                          const conv_arr_t &next_state,
                                          const conv_arr_t *output,
                                          const uint32_t msg_len,
                                          const uint8_t list_decoding_mode,
                                          const uint32_t list_size);

std::vector<std::vector<bool>> decode_post_conv_parallel_LVA(
    const std::vector<crf_mat_t> &post, const conv_arr_t &prev_state,
    const conv_arr_t &next_state, const conv_arr_t *output,
    const uint32_t msg_len, const uint8_t list_decoding_mode,
    const uint32_t list_size, const uint32_t num_thr);

std::vector<uint32_t> decode_post_vocab(
    const std::vector<crf_mat_t> &post, const uint32_t msg_len,
    const std::vector<std::vector<uint8_t>> &vocab);

std::vector<bool> viterbi_decode(const std::vector<bool> &channel_output,
                                 const conv_arr_t &prev_state,
                                 const conv_arr_t &next_state,
                                 const conv_arr_t *output,
                                 const bool must_be_perfect = false);

int main(int argc, char **argv) {
  if (argc !=4)
    throw std::runtime_error(
        "Call as ./viterbi_basecall.out decode post_file basecall_file\n\nMove array written on stdout");
  std::string mode = std::string(argv[1]);
  if (mode != "decode")
    throw std::runtime_error("invalid mode");
  std::string infile = std::string(argv[2]), outfile = std::string(argv[3]);
  /*
  if (mode == "encode") {
    // generate convolutional code matrices
    conv_arr_t prev_state, next_state, output[n_out_conv];
    generate_conv_arrays(prev_state, next_state, output);
    std::vector<bool> msg = read_bit_array(infile);
    std::vector<bool> encoded_msg = encode(msg, next_state, output);
    write_bit_array_in_bases(encoded_msg, outfile);
  }
  */
  if (mode == "decode") {
    std::vector<crf_mat_t> post = read_crf_post(infile);
    /*
    if (argc > 5) {
      if (std::string(argv[5]) != "-l") {
        std::string infile_vocab = std::string(argv[5]);
        auto vocab = read_vocab_file(infile_vocab);
        auto decoded_msg = decode_post_vocab(post, msg_len, vocab);
        write_vector(decoded_msg, outfile);
        return 0;
      }
    }
    // conv decoding
    // generate convolutional code matrices
    conv_arr_t prev_state, next_state, output[n_out_conv];
    generate_conv_arrays(prev_state, next_state, output);
    if (argc == 5) {
      std::vector<bool> decoded_msg =
          decode_post_conv(post, prev_state, next_state, output, msg_len);
      write_bit_array(decoded_msg, outfile);
    } else {
      // list decoding
      if (argc != 9 || std::string(argv[5]) != "-l")
        throw std::runtime_error("Invalid arguments.");
      uint8_t list_decoding_mode = std::stoull(std::string(argv[6]));
      uint32_t list_size = std::stoull(std::string(argv[7]));
      uint32_t num_thr = std::stoull(std::string(argv[8]));
      // do list decoding
      //
      // modes:
      //
      // 0: best sequence (not best path) list decoding (NOT WORKING SO WELL)
      //
      // 1: find list_size top paths using Parallel LVA as described in
      // https://github.com/shubhamchandak94/kBestViterbi/blob/master/kBestViterbi.py
      // or in ieeexplore.ieee.org/iel1/26/12514/00577040.pdf - if you get two
      // paths coming to same state at same time with same message, then keep
      // only one and take the max score - thus we try to get list_size unique
      // msg at each stage (if we don't do this, we observed that most of the
      // paths at the end correspond to the same msg) 
      //
      // 2: similar to above, but
      // instead of taking max score, we do logsumexp to capture the idea that
      // we want highest probability msg, not highest probability path - this
      // has less formal guarantees than 1.
      if (list_decoding_mode == 0) {
        auto decoded_msg_list = decode_post_conv_list(
            post, next_state, output, msg_len, list_decoding_mode, list_size);
      } else {
        auto decoded_msg_list = decode_post_conv_parallel_LVA(
            post, prev_state, next_state, output, msg_len, list_decoding_mode,
            list_size, num_thr);
        std::ofstream fout(outfile);
        for (auto decoded_msg : decoded_msg_list) {
          for (auto decoded_msg_bit : decoded_msg)
            fout << std::to_string(decoded_msg_bit);
          fout << "\n";
        }
        fout.close();
      }
    }*/
    // for testing
    std::vector<char> basecall = decode_post_no_conv(post);
    write_char_array(basecall, outfile);
  }
  // TODO: improve the command line options
  return 0;
}

void generate_conv_arrays(conv_arr_t &prev_state, conv_arr_t &next_state,
                          conv_arr_t *output) {
  for (uint32_t cur_state = 0; cur_state < nstate_conv; cur_state++) {
    next_state[cur_state][0] = (cur_state >> 1);
    next_state[cur_state][1] = ((cur_state | nstate_conv) >> 1);
    prev_state[cur_state][0] = (cur_state << 1) & (nstate_conv - 1);
    prev_state[cur_state][1] = ((cur_state << 1) | 1) & (nstate_conv - 1);
    output[0][cur_state][0] = __builtin_parity((cur_state & G[0]));
    output[0][cur_state][1] =
        __builtin_parity(((cur_state | nstate_conv) & G[0]));
    output[1][cur_state][0] = __builtin_parity((cur_state & G[1]));
    output[1][cur_state][1] =
        __builtin_parity(((cur_state | nstate_conv) & G[1]));
  }
  return;
}

std::vector<bool> encode(const std::vector<bool> &msg,
                         const conv_arr_t &next_state,
                         const conv_arr_t *output) {
  std::vector<bool> encoded_msg;
  uint32_t cur_state = initial_state_conv;
  for (bool msg_bit : msg) {
    encoded_msg.push_back(output[0][cur_state][msg_bit]);
    encoded_msg.push_back(output[1][cur_state][msg_bit]);
    cur_state = next_state[cur_state][msg_bit];
  }
  // add terminating bits
  for (uint8_t i = 0; i < mem_conv; i++) {
    encoded_msg.push_back(output[0][cur_state][0]);
    encoded_msg.push_back(output[1][cur_state][0]);
    cur_state = next_state[cur_state][0];
  }
  if (cur_state != 0) throw std::runtime_error("state after encoding not 0");
  return encoded_msg;
}

std::vector<bool> read_bit_array(const std::string &infile) {
  std::ifstream fin(infile);
  std::vector<bool> vec;
  char ch;
  while (fin >> std::noskipws >> ch) {
    switch (ch) {
      case '0':
        vec.push_back(0);
        break;
      case '1':
        vec.push_back(1);
        break;
      default:
        throw std::runtime_error("invalid character in input file");
    }
  }
  fin.close();
  return vec;
}

void write_bit_array(const std::vector<bool> &outvec,
                     const std::string &outfile) {
  std::ofstream fout(outfile);
  for (bool b : outvec) fout << (b ? '1' : '0');
  fout.close();
}

void write_char_array(const std::vector<char> &vec,
                      const std::string &outfile) {
  std::ofstream fout(outfile);
  for (char c : vec) fout << c;
  fout.close();
}

void write_bit_array_in_bases(const std::vector<bool> &outvec,
                              const std::string &outfile) {
  std::ofstream fout(outfile);
  uint32_t len = outvec.size();
  if (len % 2 != 0) throw std::runtime_error("length not even");
  for (uint32_t i = 0; i < len / 2; i++)
    fout << int2base[2 * outvec[2 * i] + outvec[2 * i + 1]];
  fout.close();
}

std::vector<std::vector<uint8_t>> read_vocab_file(const std::string &infile) {
  std::vector<std::vector<uint8_t>> vocab;
  std::unordered_map<char, uint8_t> char2int = {
      {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
  std::ifstream fin(infile);
  std::string line;
  while (std::getline(fin, line)) {
    std::vector<uint8_t> line_vector;
    for (char c : line) line_vector.push_back(char2int[c]);
    vocab.push_back(line_vector);
  }
  return vocab;
}

std::vector<crf_mat_t> read_crf_post(const std::string &infile) {
  std::ifstream fin(infile, std::ios::binary);
  std::vector<crf_mat_t> post;
  crf_mat_t post_mat;
  float val;
  fin.read((char *)&val, sizeof(float));
  while (!fin.eof()) {
    for (uint8_t i = 0; i < NBASE; i++) {
      for (uint8_t j = 0; j < nstate_crf; j++) {
        post_mat[i][j] = val;
        fin.read((char *)&val, sizeof(float));
      }
    }
    uint8_t i = NBASE;  // output state is one of the flop states now (for given
                        // input state, only one possible output flop state)
    for (uint8_t j = 0; j < nstate_crf; j++) {
      post_mat[i][j] = val;
      fin.read((char *)&val, sizeof(float));
    }
    post.push_back(post_mat);
  }
  return post;
}

std::vector<char> decode_post_no_conv(const std::vector<crf_mat_t> &post) {
  // just basecalling without convolutional code (for testing purposes)
  float INF = std::numeric_limits<float>::infinity();
  uint32_t nblk = post.size();
  std::vector<std::array<uint8_t, nstate_crf>> traceback(nblk);
  std::array<float, nstate_crf> curr_score, prev_score;
  curr_score.fill(0.0);

  // forward Viterbi pass
  for (uint32_t t = 0; t < nblk; t++) {
    prev_score = curr_score;
    // when to-state (st2) is flip
    for (uint8_t st2 = 0; st2 < NBASE; st2++) {
      curr_score[st2] = -INF;
      for (uint8_t st1 = 0; st1 < nstate_crf; st1++) {
        // st2 is next state, st1 is previous (previous can be flip or flop -
        // all 8 possible states)
        float score = prev_score[st1] + post[t][st2][st1];
        if (score > curr_score[st2]) {
          curr_score[st2] = score;
          traceback[t][st2] = st1;
        }
      }
    }
    // when to-state (st2) is flop (so only two possible input states)
    for (uint8_t st2 = NBASE; st2 < nstate_crf; st2++) {
      curr_score[st2] = -INF;
      uint8_t base = st2 - NBASE;  // so input state can be base or base+NBASE
      for (uint8_t ff_bit = 0; ff_bit < 2; ff_bit++) {
        uint8_t st1 = base + NBASE * ff_bit;
        // st2 is next state, st1 is previous
        float score =
            prev_score[st1] +
            post[t][NBASE]
                [st1];  // transitions to flop are stored in last row of matrix
        if (score > curr_score[st2]) {
          curr_score[st2] = score;
          traceback[t][st2] = st1;
        }
      }
    }
  }

  // traceback
  std::vector<uint8_t> path(nblk + 1);
  float score = -INF;
  for (uint8_t st = 0; st < nstate_crf; st++) {
    if (curr_score[st] > score) {
      score = curr_score[st];
      path[nblk] = st;
    }
  }
  for (uint32_t t = nblk; t > 0; t--) path[t - 1] = traceback[t - 1][path[t]];
  path.resize(nblk);
  return crfpath_to_basecall(path);
}

std::vector<char> crfpath_to_basecall(const std::vector<uint8_t> &path) {
  // collapse same states and convert to ACGT
  std::vector<char> basecall;
  for (size_t pos = 1; pos < path.size(); pos++) {
    if (path[pos] == path[pos - 1]) continue;
    std::cout << pos << "\n";
    basecall.push_back(int2base[path[pos] % NBASE]);
  }
  return basecall;
}

uint32_t get_state_idx(const uint32_t st_pos, const uint32_t st_conv,
                       const uint32_t st_crf) {
  return st_pos * nstate_conv * nstate_crf + st_conv * nstate_crf + st_crf;
}

uint32_t state_idx_to_pos(const uint32_t st) {
  return st / ((uint32_t)nstate_conv * nstate_crf);
}

uint8_t to_idx_crf_in_post(uint8_t st2_crf) {
  // return index of the st2_crf state in the post matrix, we need this because
  // we have 5 by 8 matrix and transitions to flop states are stored in the last
  // row to save space since not all transitions to the flop state are allowed
  return (st2_crf >= NBASE) ? NBASE : st2_crf;
}

std::vector<uint32_t> decode_post_vocab(
    const std::vector<crf_mat_t> &post, const uint32_t msg_len,
    const std::vector<std::vector<uint8_t>> &vocab) {
  float INF = std::numeric_limits<float>::infinity();
  uint32_t nstate_init = nstate_crf;  // initial states (same as CRF states)
  // these represent stuff before the first word occurs and can only transition
  // to the same state or to a position 0 state
  uint32_t nstate_pos = msg_len;  // number of states denoting the position
  std::vector<uint32_t> wordlen_vocab;
  for (auto word : vocab) wordlen_vocab.push_back(word.size());
  uint32_t nstate_vocab =
      std::accumulate(wordlen_vocab.begin(), wordlen_vocab.end(),
                      0);  // sum of lengths of words
  uint64_t nstate_total_64 =
      nstate_pos * nstate_vocab * 2 + nstate_init;  // 2 for flip/flop
  if (nstate_total_64 >= ((uint64_t)1 << 32))
    throw std::runtime_error("Too many states, can't fit in 32 bits");
  uint32_t nstate_total = (uint32_t)nstate_total_64;
  uint32_t nblk = post.size();
  if (post.size() < msg_len) throw std::runtime_error("Too small post matrix");

  // create unordered_map from state tuple to idx and back (does not include
  // init states, which are the first nstate_crf states)
  std::unordered_map<
      uint32_t,
      std::unordered_map<
          uint32_t,
          std::unordered_map<uint32_t, std::unordered_map<uint32_t, uint32_t>>>>
      state_tuple2idx;
  std::unordered_map<uint32_t, std::array<uint32_t, 4>> state_idx2tuple;
  uint32_t st = nstate_crf;  // first nstate_crf states are the init states
  for (uint32_t pos = 0; pos < nstate_pos; pos++) {
    for (uint32_t word_idx = 0; word_idx < wordlen_vocab.size(); word_idx++) {
      for (uint32_t pos_in_word = 0; pos_in_word < wordlen_vocab[word_idx];
           pos_in_word++) {
        for (uint32_t flip_flop_bit = 0; flip_flop_bit < 2; flip_flop_bit++) {
          state_tuple2idx[pos][word_idx][pos_in_word][flip_flop_bit] = st;
          state_idx2tuple[st] = {pos, word_idx, pos_in_word, flip_flop_bit};
          st++;
        }
      }
    }
  }

  if (st != nstate_total)
    throw std::runtime_error("Number of states don't match");

  std::vector<std::vector<uint32_t>> traceback(
      nblk, std::vector<uint32_t>(nstate_total));
  std::vector<float> curr_score(nstate_total, -INF), prev_score(nstate_total);

  // only init states allowed at the beginning, so the actual data starts at
  // first transition out of init state
  for (uint32_t st = 0; st < nstate_crf; st++) curr_score[st] = 0.0;

  // forward Viterbi pass
  for (uint32_t t = 0; t < nblk; t++) {
    prev_score = curr_score;
    // st2 is next state, st1 is previous
    // First handle st2 being one of the init states, in this case only one
    // transition allowed
    for (uint32_t st2 = 0; st2 < nstate_crf; st2++) {
      uint32_t st1 = st2;
      traceback[t][st2] = st1;
      curr_score[st2] = prev_score[st1] + post[t][to_idx_crf_in_post(st2)][st1];
    }

    // Now do all other states
    for (uint32_t st2_pos = 0; st2_pos < nstate_pos; st2_pos++) {
      for (uint32_t st2_word_idx = 0; st2_word_idx < wordlen_vocab.size();
           st2_word_idx++) {
        for (uint32_t st2_pos_in_word = 0;
             st2_pos_in_word < wordlen_vocab[st2_word_idx]; st2_pos_in_word++) {
          for (uint32_t st2_flip_flop_bit = 0; st2_flip_flop_bit < 2;
               st2_flip_flop_bit++) {
            uint32_t st2 = state_tuple2idx[st2_pos][st2_word_idx]
                                          [st2_pos_in_word][st2_flip_flop_bit];
            uint8_t st2_crf = vocab[st2_word_idx][st2_pos_in_word] +
                              st2_flip_flop_bit * NBASE;
            curr_score[st2] = -INF;

            // first, transitions from the exact same state
            uint32_t st1 = st2;
            uint8_t st1_crf = st2_crf;
            traceback[t][st2] = st1;
            curr_score[st2] =
                prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1_crf];

            // now look at transitions from init states if pos = pos_in_word = 0
            if (st2_pos == 0 && st2_pos_in_word == 0) {
              for (uint32_t st1 = 0; st1 < nstate_crf; st1++) {
                if (st1 == st2_crf)
                  continue;  // not allowed since no change in base
                if (st2_crf >= NBASE && st1 != (uint8_t)(st2_crf - NBASE))
                  continue;  // can only enter flop base from
                // same flip base
                float score =
                    prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1];
                if (score > curr_score[st2]) {
                  curr_score[st2] = score;
                  traceback[t][st2] = st1;
                }
              }
            }

            // now look at transitions from end of other words if pos != 0 and
            // pos_in_word = 0
            if (st2_pos != 0 && st2_pos_in_word == 0) {
              uint32_t st1_pos = st2_pos - 1;
              for (uint32_t st1_word_idx = 0;
                   st1_word_idx < wordlen_vocab.size(); st1_word_idx++) {
                uint32_t st1_pos_in_word = wordlen_vocab[st1_word_idx] - 1;
                for (uint32_t st1_flip_flop_bit = 0; st1_flip_flop_bit < 2;
                     st1_flip_flop_bit++) {
                  uint32_t st1 =
                      state_tuple2idx[st1_pos][st1_word_idx][st1_pos_in_word]
                                     [st1_flip_flop_bit];
                  uint8_t st1_crf = vocab[st1_word_idx][st1_pos_in_word] +
                                    st1_flip_flop_bit * NBASE;
                  if (st1_crf == st2_crf)
                    continue;  // not allowed since no change in base
                  if (st2_crf >= NBASE && st1_crf != st2_crf - NBASE)
                    continue;  // can only enter flop base from
                  // same flip base
                  float score = prev_score[st1] +
                                post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
                  if (score > curr_score[st2]) {
                    curr_score[st2] = score;
                    traceback[t][st2] = st1;
                  }
                }
              }
            }

            // finally look at transitions from same word if pos_in_word != 0
            if (st2_pos_in_word != 0) {
              uint32_t st1_pos = st2_pos;
              uint32_t st1_word_idx = st2_word_idx;
              uint32_t st1_pos_in_word = st2_pos_in_word - 1;
              for (uint32_t st1_flip_flop_bit = 0; st1_flip_flop_bit < 2;
                   st1_flip_flop_bit++) {
                uint32_t st1 =
                    state_tuple2idx[st1_pos][st1_word_idx][st1_pos_in_word]
                                   [st1_flip_flop_bit];
                uint8_t st1_crf = vocab[st1_word_idx][st1_pos_in_word] +
                                  st1_flip_flop_bit * NBASE;
                if (st1_crf == st2_crf)
                  continue;  // not allowed since no change in base
                if (st2_crf >= NBASE && st1_crf != st2_crf - NBASE)
                  continue;  // can only enter flop base from
                // same flip base
                float score = prev_score[st1] +
                              post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
                if (score > curr_score[st2]) {
                  curr_score[st2] = score;
                  traceback[t][st2] = st1;
                }
              }
            }
          }
        }
      }
    }
  }
  // traceback
  std::vector<uint32_t> path(nblk + 1);
  float score = -INF;
  uint32_t st_pos = nstate_pos - 1;  // pos is at end for last state
  for (uint32_t st_word_idx = 0; st_word_idx < wordlen_vocab.size();
       st_word_idx++) {
    uint32_t st_pos_in_word =
        wordlen_vocab[st_word_idx] - 1;  // last state at end of word
    for (uint32_t st_flip_flop_bit = 0; st_flip_flop_bit < 2;
         st_flip_flop_bit++) {
      uint32_t st = state_tuple2idx[st_pos][st_word_idx][st_pos_in_word]
                                   [st_flip_flop_bit];
      if (curr_score[st] > score) {
        score = curr_score[st];
        path[nblk] = st;
      }
    }
  }

  for (uint32_t t = nblk; t > 0; t--) path[t - 1] = traceback[t - 1][path[t]];
  // decode the msg from the state sequence
  std::vector<uint32_t> decoded_msg;
  int64_t cur_pos = -1;
  for (uint32_t t = 0; t <= nblk; t++) {
    // look for non-init states where pos increases
    uint32_t st = path[t];
    if (st < nstate_crf) continue;  // init state, skip
    auto st_tuple = state_idx2tuple[st];
    if (st_tuple[0] > cur_pos) {
      if (st_tuple[0] != cur_pos + 1)
        throw std::runtime_error("pos increase not 1");
      if (st_tuple[2] != 0)
        throw std::runtime_error("pos_in_word at transition not 0");
      cur_pos = st_tuple[0];
      decoded_msg.push_back(st_tuple[1]);
    }
  }
  if (decoded_msg.size() != msg_len)
    throw std::runtime_error("Decoded message length does not match msg_len");
  return decoded_msg;
}

std::vector<path_t> decode_post_conv_list(const std::vector<crf_mat_t> &post,
                                          const conv_arr_t &next_state,
                                          const conv_arr_t *output,
                                          const uint32_t msg_len,
                                          const uint8_t list_decoding_mode,
                                          const uint32_t list_size) {
  if (post.size() < msg_len + mem_conv)
    throw std::runtime_error("Too small post matrix");
  uint32_t nblk = post.size();
  if (list_decoding_mode != 0)
    throw std::runtime_error(
        "Invalid list decoding mode");  // currently only mode 0 supported
  std::vector<path_t> prev_path_list, curr_path_list;
  // initialize path list
  // states with msg = '', pos = 0, st_conv = initial_state_conv, all st_crf's,
  // score = 0
  for (uint8_t st_crf = 0; st_crf < nstate_crf; st_crf++)
    curr_path_list.push_back({std::string(""), 0, initial_state_conv, st_crf});

  for (uint32_t t = 0; t < nblk; t++) {
    prev_path_list = curr_path_list;
    curr_path_list.clear();
    // create maps from msg to position in curr_list
    // (need two maps since msg can correspond to 2 states in some cases
    // basically the difference between A+A-A+ and A-A+A- which have
    // same message but end up at different CRF states
    // These maps are used for combining paths (where relevant)
    std::unordered_map<std::string, uint32_t>
        msg2listpos[2];  // 0 - flip, 1 - flop
    std::vector<path_t>
        new_path_list;  // stores new paths with possible duplicates (resolved
                        // when inserting to curr_path_list)
    path_t new_path;
    // now go through the paths one by one, extending them in all possible ways
    for (auto prev_path : prev_path_list) {
      // first do transition to same state
      new_path.msg = prev_path.msg;
      new_path.st_pos = prev_path.st_pos;
      new_path.st_conv = prev_path.st_conv;
      new_path.st_crf = prev_path.st_crf;
      new_path.score =
          prev_path.score +
          post[t][to_idx_crf_in_post(new_path.st_crf)][prev_path.st_crf];
      new_path_list.push_back(new_path);

      if (prev_path.st_pos ==
          msg_len +
              mem_conv)  // at final stage, so no more next states to go to.
        continue;

      // now transition to next state when next input is 0/1
      for (uint8_t conv_bit = 0; conv_bit < 2; conv_bit++) {
        // first check if that's not allowed due to sync constraints
        // note that st_pos is sort of 1 indexed and so prev_path.st_pos
        // actually represents the position of the upcoming  conv_bit
        if (prev_path.st_pos < msg_len &&
            prev_path.st_pos % sync_marker_period < sync_marker_length)
          if (conv_bit != sync_marker[prev_path.st_pos % sync_marker_period])
            continue;

        // if we are at msg_len, bit must be 0
        if (prev_path.st_pos >= msg_len && conv_bit == 1) continue;

        new_path.msg = prev_path.msg + std::to_string(conv_bit);
        new_path.st_pos = prev_path.st_pos + 1;
        new_path.st_conv = next_state[prev_path.st_conv][conv_bit];
        uint8_t prev_st_crf_base = prev_path.st_crf % NBASE;
        uint8_t new_st_crf_base = 2 * output[0][prev_path.st_conv][conv_bit] +
                                  output[1][prev_path.st_conv][conv_bit];
        // next crf_state depends on whether the prev and new bases are same
        if (new_st_crf_base != prev_st_crf_base)  //  going to flip state
          new_path.st_crf = new_st_crf_base;
        else  // flip or flop according to whether current is flop or flip
          new_path.st_crf = (prev_path.st_crf < NBASE)
                                ? (prev_path.st_crf + NBASE)
                                : (prev_path.st_crf - NBASE);
        new_path.score =
            prev_path.score +
            post[t][to_idx_crf_in_post(new_path.st_crf)][prev_path.st_crf] + 20;
        new_path_list.push_back(new_path);
      }
    }
    // now put the paths into curr_path_list merging when relevant
    for (auto new_path : new_path_list) {
      // check if already exists, if so combine
      if (msg2listpos[(new_path.st_crf < NBASE)].count(new_path.msg) == 1) {
        uint32_t pos_in_curr =
            msg2listpos[(new_path.st_crf < NBASE)][new_path.msg];
        curr_path_list[pos_in_curr]
            .score =  // std::max(new_path.score,curr_path_list[pos_in_curr].score);
            logsumexpf(new_path.score, curr_path_list[pos_in_curr].score);
      } else {
        msg2listpos[(new_path.st_crf < NBASE)][new_path.msg] =
            curr_path_list.size();
        curr_path_list.push_back(new_path);
      }
    }

    // finally take only top list_size paths in the list
    std::nth_element(curr_path_list.begin(), curr_path_list.begin() + list_size,
                     curr_path_list.end(),
                     [](const path_t &p1, const path_t &p2) -> bool {
                       return p1.score > p2.score;
                     });
    //        std::sort(curr_path_list.begin(),curr_path_list.end(),[](const
    //        path_t &p1, const path_t &p2) -> bool {
    //         return p1.score > p2.score;
    //        });
    if (curr_path_list.size() > list_size) curr_path_list.resize(list_size);

    // FOR DEBUGGING: PRINT THE STATES IN THE LIST
    /*
    std::cout << "t: " << t << "\n";
    for (auto curr_path : curr_path_list) {
      std::cout << curr_path.msg << "\n"
                << curr_path.st_pos << "\n"
                << curr_path.score << "\n";
      std::cout << "\n";
    }
    */
  }
  return curr_path_list;
}

std::vector<bool> decode_post_conv(const std::vector<crf_mat_t> &post,
                                   const conv_arr_t &prev_state,
                                   const conv_arr_t &next_state,
                                   const conv_arr_t *output,
                                   const uint32_t msg_len) {
  float INF = std::numeric_limits<float>::infinity();
  uint32_t nstate_pos =
      msg_len + mem_conv +
      1;  // number of states denoting the position in convolution trellis
  // 0 to msg_len + mem_conv
  uint64_t nstate_total_64 = nstate_pos * nstate_crf * nstate_conv;
  if (nstate_total_64 >= ((uint64_t)1 << 32))
    throw std::runtime_error("Too many states, can't fit in 32 bits");
  uint32_t nstate_total = (uint32_t)nstate_total_64;
  uint32_t nblk = post.size();
  if (post.size() < msg_len + mem_conv)
    throw std::runtime_error("Too small post matrix");
  std::vector<std::vector<uint32_t>> traceback(
      nblk, std::vector<uint32_t>(nstate_total));
  std::vector<float> curr_score(nstate_total, -INF), prev_score(nstate_total);

  for (uint8_t st_crf = 0; st_crf < nstate_crf; st_crf++)
    curr_score[get_state_idx(0, initial_state_conv, st_crf)] =
        0.0;  // only valid initial state is pos 0, conv code at
              // initial_state_conv. crf state can be anything since we will
              // later ignore everything before the first transition

  // forward Viterbi pass
  for (uint32_t t = 0; t < nblk; t++) {
    prev_score = curr_score;
    // st2 is next state, st1 is previous
    for (uint32_t st2_pos = 0; st2_pos < nstate_pos; st2_pos++) {
      for (uint32_t st2_conv = 0; st2_conv < nstate_conv; st2_conv++) {
        for (uint8_t st2_crf = 0; st2_crf < nstate_crf; st2_crf++) {
          uint32_t st2 = get_state_idx(st2_pos, st2_conv, st2_crf);
          curr_score[st2] = -INF;
          // now we have two possibilities, st2_crf can be flip or flop state.
          // if flip, then we go through all input states, and except for the
          // case when st1_crf == st2_crf, previous state has one less pos and
          // may or may have a valid previous conv state depending on the base
          // in question if st2_crf is flop, then there are only two possible
          // input states the flip and flop state for the same base. For the
          // flop state, everything else stays the same. For the flip input
          // state, it has one less pos (and so on).

          for (uint8_t st1_crf = 0; st1_crf < nstate_crf; st1_crf++) {
            if (st2_crf >= NBASE &&
                !((st1_crf == st2_crf) || st1_crf == st2_crf - NBASE))
              continue;  // unallowed transition
            if (st2_crf == st1_crf) {
              // at same base
              uint32_t st1 = st2;
              float score = prev_score[st1] +
                            post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
              if (score > curr_score[st2]) {
                curr_score[st2] = score;
                traceback[t][st2] = st1;
              }
            } else {
              // new base, new position
              // look at two possible previous states of convolutional code and
              // see if the output for the transition matches the base st2_crf
              bool curr_conv_bit = (st2_conv >> (mem_conv - 1));
              uint8_t st2_crf_base = st2_crf % NBASE;
              if (st2_pos == 0)
                continue;  // can't go to new state while still being at
                           // position 0
              uint32_t st1_pos = st2_pos - 1;

              // sync_markers
              if (st1_pos < msg_len &&
                  (st1_pos % sync_marker_period < sync_marker_length))
                if (curr_conv_bit != sync_marker[st1_pos % sync_marker_period])
                  continue;  // invalid transition
              for (uint8_t conv_bit = 0; conv_bit < 2; conv_bit++) {
                uint32_t st1_conv = prev_state[st2_conv][conv_bit];
                if (2 * output[0][st1_conv][curr_conv_bit] +
                        output[1][st1_conv][curr_conv_bit] ==
                    st2_crf_base) {
                  uint32_t st1 = get_state_idx(st1_pos, st1_conv, st1_crf);
                  float score = prev_score[st1] +
                                post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
                  if (score > curr_score[st2]) {
                    curr_score[st2] = score;
                    traceback[t][st2] = st1;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // traceback
  std::vector<uint32_t> path(nblk + 1);
  std::vector<uint8_t> crfpath(nblk + 1);
  float score = -INF;
  uint32_t st_pos = msg_len + mem_conv, st_conv = 0;  // last state
  for (uint8_t st_crf = 0; st_crf < nstate_crf; st_crf++) {
    uint32_t st = get_state_idx(st_pos, st_conv, st_crf);
    if (curr_score[st] > score) {
      score = curr_score[st];
      path[nblk] = st;
    }
  }
  std::cout << "Best path score: " << score << "\n";
  for (uint32_t t = nblk; t > 0; t--) path[t - 1] = traceback[t - 1][path[t]];
  for (uint32_t t = 0; t < nblk + 1; t++) crfpath[t] = path[t] % nstate_crf;
  std::vector<char> basecall = crfpath_to_basecall(crfpath);
  if (basecall.size() != msg_len + mem_conv)
    throw std::runtime_error("incorrect decoded length");
  // convert basecall to bool vector and then decode using viterbi decode
  // function
  std::vector<bool> channel_output;
  for (char c : basecall) {
    switch (c) {
      case 'A':
        channel_output.push_back(0);
        channel_output.push_back(0);
        break;
      case 'C':
        channel_output.push_back(0);
        channel_output.push_back(1);
        break;
      case 'G':
        channel_output.push_back(1);
        channel_output.push_back(0);
        break;
      case 'T':
        channel_output.push_back(1);
        channel_output.push_back(1);
        break;
      default:
        throw std::runtime_error("unexpected character in basecall");
    }
  }
  return viterbi_decode(channel_output, prev_state, next_state, output, true);
}

float combine_scores(const float &score_1, const float &score_2, const uint8_t &list_decoding_mode) {
    if (list_decoding_mode == 1) 
        return std::max(score_1,score_2);
    else
        return logsumexpf(score_1,score_2);
}

std::vector<std::vector<bool>> decode_post_conv_parallel_LVA(
    const std::vector<crf_mat_t> &post, const conv_arr_t &prev_state,
    const conv_arr_t &next_state, const conv_arr_t *output,
    const uint32_t msg_len, const uint8_t list_decoding_mode,
    const uint32_t list_size, const uint32_t num_thr) {

  omp_set_num_threads(num_thr);
  float INF = std::numeric_limits<float>::infinity();
  uint32_t nstate_pos =
      msg_len + mem_conv +
      1;  // number of states denoting the position in convolution trellis
  // 0 to msg_len + mem_conv
  uint64_t nstate_total_64 = nstate_pos * nstate_crf * nstate_conv;
  if (nstate_total_64 >= ((uint64_t)1 << 32))
    throw std::runtime_error("Too many states, can't fit in 32 bits");
  uint32_t nstate_total = (uint32_t)nstate_total_64;
  uint32_t nblk = post.size();
  if (post.size() < msg_len + mem_conv)
    throw std::runtime_error("Too small post matrix");

  // instead of traceback, store the msg till now as a bitset
  if (msg_len > BITSET_SIZE)
    throw std::runtime_error("msg_len can't be above BITSET_SIZE");

  bitset_t **curr_best_msg = new bitset_t *[list_size];
  bitset_t **prev_best_msg = new bitset_t *[list_size];
  float **curr_score = new float *[list_size];
  float **prev_score = new float *[list_size];
  for (uint32_t i = 0; i < list_size; i++) {
    curr_best_msg[i] = new bitset_t[nstate_total];
    prev_best_msg[i] = new bitset_t[nstate_total];
    curr_score[i] = new float[nstate_total];
    prev_score[i] = new float[nstate_total];
    for (uint32_t j = 0; j < nstate_total; j++) {
      curr_score[i][j] = -INF;
      prev_score[i][j] = -INF;
    }
  }

  // as opposed to best path Viterbi, now we first store all
  // incoming paths in an array
  LVA_path_t **LVA_path_list = new LVA_path_t * [num_thr];
  for (uint32_t i = 0; i < num_thr; i++) {
      LVA_path_list[i] = new LVA_path_t[15*list_size];
      // at most 15 incoming transitions for a state
  }

  for (uint8_t st_crf = 0; st_crf < nstate_crf; st_crf++) {
    curr_score[0][get_state_idx(0, initial_state_conv, st_crf)] =
        0.0;  // only valid initial state is pos 0, conv code at
              // initial_state_conv. crf state can be anything since we will
              // later ignore everything before the first transition
              // only populate one position in list (0), rest -INF
  }

  // forward Viterbi pass
  for (uint32_t t = 0; t < nblk; t++) {
    // swap prev and curr arrays
    std::swap(curr_score, prev_score);
    std::swap(curr_best_msg, prev_best_msg);

// st2 is next state, st1 is previous
    for (uint32_t st2_pos =
             std::max((int64_t)nstate_pos - 2 - (nblk - 1 - t), (int64_t)0);
         st2_pos < std::min(t + 2, nstate_pos); st2_pos++) {
      // only allow pos which can have non -INF scores or will lead to useful
      // final states initially large pos is not allowed, and at the end small
      // pos not allowed (since those can't lead to correct st2_pos at the end).
#pragma omp parallel
#pragma omp for
      for (uint32_t st2_conv = 0; st2_conv < nstate_conv; st2_conv++) {
        uint32_t tid = omp_get_thread_num();
        for (uint8_t st2_crf = 0; st2_crf < nstate_crf; st2_crf++) {
          bool curr_conv_bit = (st2_conv >> (mem_conv - 1));
          if (st2_pos > msg_len && curr_conv_bit == 1)
            continue;  // invalid state since padding is all 0's
          uint32_t st2 = get_state_idx(st2_pos, st2_conv, st2_crf);
          uint32_t pos_in_Lpl = 0;
          // we will later on deduplicate LVA_path_list - the only way we
          // can have a duplicate is when one of the stay transition msg
          // matches one of the non-stay transitions or between non-stay
          // transitions coming from same base (one flip, one flop)
          // note that stuff coming from same state can't be duplicates since
          // we also deduplicated at previous step
          uint32_t startpos_stay = 0, endpos_stay, startpos_flip[NBASE], endpos_flip[NBASE];
          
          // now we have two possibilities, st2_crf can be flip or flop state.
          // if flip, then we go through all input states, and except for the
          // case when st1_crf == st2_crf, previous state has one less pos and
          // may or may have a valid previous conv state depending on the base
          // in question if st2_crf is flop, then there are only two possible
          // input states the flip and flop state for the same base. For the
          // flop state, everything else stays the same. For the flip input
          // state, it has one less pos (and so on).
        
          // first do stay
          endpos_stay = 0;
          uint8_t st1_crf = st2_crf;
          uint32_t st1 = st2;
          for (uint32_t list_pos = 0; list_pos < list_size; list_pos++) {
            if (prev_score[list_pos][st1] == -INF) break;
            float score = prev_score[list_pos][st1] +
                      post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
            LVA_path_list[tid][pos_in_Lpl++] = {prev_best_msg[list_pos][st1], score};
            endpos_stay++;
          }
          for (uint8_t st1_crf = 0; st1_crf < nstate_crf; st1_crf++) {
            if (st2_crf >= NBASE &&
                !((st1_crf == st2_crf) || st1_crf == st2_crf - NBASE))
              continue;  // unallowed transition
            if (st2_crf == st1_crf) {
              // at same base, already done above
              continue;
            } else {
              // new base, new position
              // look at two possible previous states of convolutional code and
              // see if the output for the transition matches the base st2_crf

              uint8_t st2_crf_base = st2_crf % NBASE;
              if (st2_pos == 0)
                continue;  // can't go to new state while still being at
                           // position 0
              uint32_t st1_pos = st2_pos - 1;

              // sync_markers
              if (st1_pos < msg_len &&
                  (st1_pos % sync_marker_period < sync_marker_length))
                if (curr_conv_bit != sync_marker[st1_pos % sync_marker_period])
                  continue;  // invalid transition
              if (st1_crf < NBASE)
                  endpos_flip[st1_crf] = startpos_flip[st1_crf] = pos_in_Lpl;
              for (uint8_t conv_bit = 0; conv_bit < 2; conv_bit++) {
                uint32_t st1_conv = prev_state[st2_conv][conv_bit];
                if (2 * output[0][st1_conv][curr_conv_bit] +
                        output[1][st1_conv][curr_conv_bit] ==
                    st2_crf_base) {
                  uint32_t st1 = get_state_idx(st1_pos, st1_conv, st1_crf);
                  for (uint32_t list_pos = 0; list_pos < list_size;
                       list_pos++) {
                    if (prev_score[list_pos][st1] == -INF) break;
                    float score = prev_score[list_pos][st1] +
                                  post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
                    bitset_t next_msg;
                    if (st2_pos > msg_len) {
                      next_msg = prev_best_msg[list_pos]
                                              [st1];  // no need to push in 0 in
                                                      // this case since padding
                                                      // is known to be 0
                    } else {
                      next_msg = (prev_best_msg[list_pos][st1] << 1);
                      next_msg |= curr_conv_bit;
                    }
                    LVA_path_t temp_path = {next_msg,0};
                    bool foundmatch = false;
                    auto findpos = std::find(LVA_path_list[tid]+startpos_stay,LVA_path_list[tid]+endpos_stay,temp_path);
                    if (findpos != LVA_path_list[tid]+endpos_stay) {
                        foundmatch = true;
                    } else if (st1_crf >= NBASE) {
                        // not found and st1 is flop (so search in the flip)
                        findpos = std::find(LVA_path_list[tid]+startpos_flip[st1_crf-NBASE],LVA_path_list[tid]+endpos_flip[st1_crf-NBASE],temp_path);
                        if (findpos != LVA_path_list[tid]+endpos_flip[st1_crf-NBASE])
                            foundmatch = true;
                    }

                    if (foundmatch == true) {
                        // found match
                        findpos->score = combine_scores(score,findpos->score,list_decoding_mode);
                    } else {
                        LVA_path_list[tid][pos_in_Lpl++] = {next_msg, score};
                        if (st1_crf < NBASE)
                            endpos_flip[st1_crf]++; // flip base
                    }
                  }
                }
              }
            }
          }

          // pick top list_size elements from LVA_path_list
          if (pos_in_Lpl > list_size)
            std::nth_element(
                LVA_path_list[tid],
                LVA_path_list[tid] + list_size,
                LVA_path_list[tid] + pos_in_Lpl,
                [](const LVA_path_t &p1, const LVA_path_t &p2) -> bool {
                  return p1.score > p2.score;
                });

          // update curr_best_msg and curr_score (fill rest with -INF)
          for (uint32_t list_pos = 0; list_pos < list_size; list_pos++) {
            if (list_pos < pos_in_Lpl) {
              curr_score[list_pos][st2] = LVA_path_list[tid][list_pos].score;
              curr_best_msg[list_pos][st2] = LVA_path_list[tid][list_pos].msg;
            } else {
              curr_score[list_pos][st2] = -INF;
            }
          }
        } 
      }
    }
  }
  std::vector<LVA_path_t> LVA_path_list_final;
  uint32_t st_pos = msg_len + mem_conv, st_conv = 0;  // last state
  for (uint8_t st_crf = 0; st_crf < nstate_crf; st_crf++) {
    uint32_t st = get_state_idx(st_pos, st_conv, st_crf);
    for (uint32_t list_pos = 0; list_pos < list_size; list_pos++) {
      if (curr_score[list_pos][st] != -INF)
        LVA_path_list_final.push_back(
            {curr_best_msg[list_pos][st], curr_score[list_pos][st]});
    }
  }

  // sort
  std::sort(LVA_path_list_final.begin(), LVA_path_list_final.end(),
            [](const LVA_path_t &p1, const LVA_path_t &p2) -> bool {
              return p1.score > p2.score;
            });

  std::vector<std::vector<bool>> decoded_msg_list;

  // now convert bitset to bool vectors
  for (auto LVA_path : LVA_path_list_final) {
    std::vector<bool> decoded_msg(msg_len);
    for (uint8_t i = 0; i < msg_len; i++)
      decoded_msg[i] =
          LVA_path
              .msg[msg_len - 1 - i];  // due to way bitset is stored in reverse
    decoded_msg_list.push_back(decoded_msg);
    // FOR DEBUGGING
/*
    std::cout << "score: " << LVA_path.score << "\n";
    for (auto b : decoded_msg_list.back()) std::cout << b;
    std::cout << "\n\n";
*/    
  }
//  std::cout << "Final list size: " << decoded_msg_list.size() << "\n";

  for (uint8_t i = 0; i < list_size; i++) {
    delete[] curr_best_msg[i];
    delete[] prev_best_msg[i];
    delete[] curr_score[i];
    delete[] prev_score[i];
  }
  delete[] curr_best_msg;
  delete[] prev_best_msg;
  delete[] curr_score;
  delete[] prev_score;
  for (uint32_t i = 0; i < num_thr; i++)
      delete[] LVA_path_list[i];
  delete[] LVA_path_list;
  return decoded_msg_list;
}

std::vector<bool> viterbi_decode(const std::vector<bool> &channel_output,
                                 const conv_arr_t &prev_state,
                                 const conv_arr_t &next_state,
                                 const conv_arr_t *output,
                                 const bool must_be_perfect) {
  // must_be_perfect flag for cases when we expect 0 errors (e.g., when called
  // from decode_post_conv)
  float INF = std::numeric_limits<float>::infinity();
  uint32_t out_size = channel_output.size();
  if (out_size % n_out_conv != 0)
    throw std::runtime_error("length not multiple of n_out_conv");
  uint32_t in_size = out_size / n_out_conv;
  if (in_size < (uint32_t)mem_conv)
    throw std::runtime_error("too small channel output");
  std::vector<std::array<uint32_t, nstate_conv>> traceback(in_size);
  std::vector<float> curr_score(nstate_conv, -INF), prev_score(nstate_conv);
  curr_score[initial_state_conv] = 0.0;
  for (uint32_t t = 0; t < in_size; t++) {
    prev_score = curr_score;
    for (uint32_t st2 = 0; st2 < nstate_conv; st2++) {
      // st2 = next state
      uint32_t st1 = prev_state[st2][0];
      bool curr_bit = (st2 >> (mem_conv - 1));
      // sync_markers
      if (t < in_size - mem_conv &&
          (t % sync_marker_period < sync_marker_length))
        if (curr_bit != sync_marker[t % sync_marker_period])
          continue;  // invalid transition
      curr_score[st2] =
          prev_score[st1] -
          (float)(channel_output[2 * t] != output[0][st1][curr_bit]) -
          (float)(channel_output[2 * t + 1] != output[1][st1][curr_bit]);
      traceback[t][st2] = st1;
      st1 = prev_state[st2][1];
      float score =
          prev_score[st1] -
          (float)(channel_output[2 * t] != output[0][st1][curr_bit]) -
          (float)(channel_output[2 * t + 1] != output[1][st1][curr_bit]);
      if (score > curr_score[st2]) {
        curr_score[st2] = score;
        traceback[t][st2] = st1;
      }
    }
  }
  if (must_be_perfect && (curr_score[0] != 0.0))
    throw std::runtime_error(
        "found errors in output sequence even though none were expected");
  std::vector<bool> decoded_msg(in_size);
  uint32_t cur_state = 0;  // we already know the last state is 0
  decoded_msg[in_size - 1] = (cur_state >> (mem_conv - 1));
  for (uint32_t t = in_size - 1; t > 0; t--) {
    cur_state = traceback[t][cur_state];
    decoded_msg[t - 1] = (cur_state >> (mem_conv - 1));
  }
  decoded_msg.resize(in_size - mem_conv);
  return decoded_msg;
}
