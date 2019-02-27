#include <array>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <numeric>

const uint8_t NBASE = 4;
const char int2base[NBASE] = {'A', 'C', 'G', 'T'};
const bool base2bit[NBASE][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
const uint8_t nstate_crf = 8;
// states are A+,C+,G+,T+,A-,C-,G-,T- (+ = flip, - = flop), you can enter from a
// different base only into a flip base.

typedef std::array<std::array<float, nstate_crf>, NBASE + 1> crf_mat_t;
// for each input state, we have five (4+1) possible output states: e.g., A+ ->
// A+,C+,G+,T+,A-; A- -> A+,C+,G+,T+,A-; A-; etc.

// convolutional code related parameters
const uint8_t mem_conv =
    // 6; // mem 6 from CCSDS
//        8;  // mem 8 from GL paper
    11;  // mem 11 from GL paper
// 14; // mem 14 from GL paper
const uint32_t nstate_conv = 1 << mem_conv;
const uint8_t n_out_conv = 2;
typedef std::array<std::array<uint32_t, 2>, nstate_conv> conv_arr_t;
const uint32_t G[n_out_conv] =  // octal
                                // {0171, 0133};  // mem 6 from CCSDS
//        {0515, 0677};  // mem 8 from GL paper
    {05537, 06131};  // mem 11 from GL paper
// {075063, 056711}; // mem 14 from GL paper
const uint32_t initial_state_conv =  // binary
    //    0;                               // 0 initial state
//     0b100101; // mem 6
//     0b10010110; // mem 8
   0b10010110001;  // mem 11
// 0b10010110001101; // mem 14

// when using sync_markers, initial_state = 0 should work just fine
const uint32_t sync_marker_length = //0;
// 1
// 2;
3;
const char sync_marker[sync_marker_length] = 
// {};
// {1};
// {1,0};
{1, 1, 0};
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

std::vector<uint32_t> decode_post_vocab(const std::vector<crf_mat_t> &post,
                                   const uint32_t msg_len,
                                   const std::vector<std::vector<uint8_t>> &vocab);

std::vector<bool> viterbi_decode(const std::vector<bool> &channel_output,
                                 const conv_arr_t &prev_state,
                                 const conv_arr_t &next_state,
                                 const conv_arr_t *output,
                                 const bool must_be_perfect = false);

int main(int argc, char **argv) {

  if (argc < 4)
    throw std::runtime_error(
        "not enough arguments. Call as ./a.out [encode/decode] infile outfile "
        "[msg_len_for_decode] [infile_for_vocab]");
  std::string mode = std::string(argv[1]);
  if (mode != "encode" && mode != "decode")
    throw std::runtime_error("invalid mode");
  std::string infile = std::string(argv[2]), outfile = std::string(argv[3]);
  if (mode == "encode") {
    // generate convolutional code matrices
    conv_arr_t prev_state, next_state, output[n_out_conv];
    generate_conv_arrays(prev_state, next_state, output);
    std::vector<bool> msg = read_bit_array(infile);
    std::vector<bool> encoded_msg = encode(msg, next_state, output);
    write_bit_array_in_bases(encoded_msg, outfile);
  }
  if (mode == "decode") {
    if (argc < 5)
      throw std::runtime_error(
          "not enough arguments. Call as ./a.out [encode/decode] infile "
          "outfile [msg_len_for_decode] [infile_for_vocab]");
    uint32_t msg_len = std::stoull(std::string(argv[4]));
    std::vector<crf_mat_t> post = read_crf_post(infile);
    if (argc == 5) {
        // conv decoding
        // generate convolutional code matrices
        conv_arr_t prev_state, next_state, output[n_out_conv];
        generate_conv_arrays(prev_state, next_state, output);
        std::vector<bool> decoded_msg =
            decode_post_conv(post, prev_state, next_state, output, msg_len);
        write_bit_array(decoded_msg, outfile);
    } else {
        std::string infile_vocab = std::string(argv[5]);
        auto vocab = read_vocab_file(infile_vocab);
        auto decoded_msg = decode_post_vocab(post, msg_len, vocab);
        write_vector(decoded_msg, outfile);
    }
    // for testing
    //        std::vector<char> basecall = decode_post_no_conv(post);
    //        write_char_array(basecall, outfile);
  }
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
    std::unordered_map<char,uint8_t> char2int = {{'A',0},{'C',1},{'G',2},{'T',3}};
    std::ifstream fin(infile);
    std::string line;
    while(std::getline(fin,line)) {
        std::vector<uint8_t> line_vector;
        for (char c: line)
            line_vector.push_back(char2int[c]);
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
    // return index of the st2_crf state in the post matrix, we need this because we 
    // have 5 by 8 matrix and transitions to flop states are stored in the last row to save space
    // since not all transitions to the flop state are allowed
    return (st2_crf >= NBASE) ? NBASE : st2_crf;
} 

std::vector<uint32_t> decode_post_vocab(const std::vector<crf_mat_t> &post,
                                   const uint32_t msg_len,
                                   const std::vector<std::vector<uint8_t>> &vocab) {
  float INF = std::numeric_limits<float>::infinity();
  uint32_t nstate_init = nstate_crf; // initial states (same as CRF states)
  // these represent stuff before the first word occurs and can only transition to the 
  // same state or to a position 0 state
  uint32_t nstate_pos =
      msg_len;  // number of states denoting the position
  std::vector<uint32_t> wordlen_vocab;
  for (auto word: vocab)
      wordlen_vocab.push_back(word.size());
  uint32_t nstate_vocab = std::accumulate(wordlen_vocab.begin(),wordlen_vocab.end(),0);// sum of lengths of words
  uint64_t nstate_total_64 = nstate_pos * nstate_vocab * 2 + nstate_init; // 2 for flip/flop
  if (nstate_total_64 >= ((uint64_t)1 << 32))
    throw std::runtime_error("Too many states, can't fit in 32 bits");
  uint32_t nstate_total = (uint32_t)nstate_total_64;
  uint32_t nblk = post.size();
  if (post.size() < msg_len)
    throw std::runtime_error("Too small post matrix");

  // create unordered_map from state tuple to idx and back (does not include init states, which are the first nstate_crf states)
  std::unordered_map<uint32_t,std::unordered_map<uint32_t,std::unordered_map<uint32_t,std::unordered_map<uint32_t,uint32_t>>>> state_tuple2idx;
  std::unordered_map<uint32_t,std::array<uint32_t,4>> state_idx2tuple;
  uint32_t st = nstate_crf; // first nstate_crf states are the init states
  for (uint32_t pos = 0; pos < nstate_pos; pos++) {
      for (uint32_t word_idx = 0; word_idx < wordlen_vocab.size(); word_idx++) {
          for (uint32_t pos_in_word = 0; pos_in_word < wordlen_vocab[word_idx]; pos_in_word++) {
              for (uint32_t flip_flop_bit = 0; flip_flop_bit < 2; flip_flop_bit++) {
                  state_tuple2idx[pos][word_idx][pos_in_word][flip_flop_bit] = st;
                  state_idx2tuple[st] = {pos,word_idx,pos_in_word,flip_flop_bit};
                  st++;
              }
          }
      }
  }

  if (st != nstate_total) throw std::runtime_error("Number of states don't match");

  std::vector<std::vector<uint32_t>> traceback(
      nblk, std::vector<uint32_t>(nstate_total));
  std::vector<float> curr_score(nstate_total, -INF), prev_score(nstate_total);

  // only init states allowed at the beginning, so the actual data starts at first transition out of init state
  for (uint32_t st = 0; st < nstate_crf; st++)
      curr_score[st] = 0.0;

  // forward Viterbi pass
  for (uint32_t t = 0; t < nblk; t++) {
    prev_score = curr_score;
    // st2 is next state, st1 is previous
    // First handle st2 being one of the init states, in this case only one transition allowed
    for (uint32_t st2 = 0; st2 < nstate_crf; st2++) {
        uint32_t st1 = st2;
        traceback[t][st2] = st1;
        curr_score[st2] = prev_score[st1] + post[t][to_idx_crf_in_post(st2)][st1];
    }
    
    // Now do all other states
    for (uint32_t st2_pos = 0; st2_pos < nstate_pos; st2_pos++) {
        for (uint32_t st2_word_idx = 0; st2_word_idx < wordlen_vocab.size(); st2_word_idx++) {
            for (uint32_t st2_pos_in_word = 0; st2_pos_in_word < wordlen_vocab[st2_word_idx]; st2_pos_in_word++) {
                for (uint32_t st2_flip_flop_bit = 0; st2_flip_flop_bit < 2; st2_flip_flop_bit++) {
                    uint32_t st2 = state_tuple2idx[st2_pos][st2_word_idx][st2_pos_in_word][st2_flip_flop_bit];
                    uint8_t st2_crf = vocab[st2_word_idx][st2_pos_in_word] + st2_flip_flop_bit*NBASE;
                    curr_score[st2] = -INF;

                    // first, transitions from the exact same state
                    uint32_t st1 = st2;
                    uint8_t st1_crf = st2_crf;
                    traceback[t][st2] = st1;
                    curr_score[st2] = prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
                    
                    // now look at transitions from init states if pos = pos_in_word = 0
                    if (st2_pos == 0 && st2_pos_in_word == 0) {
                        for (uint32_t st1 = 0; st1 < nstate_crf; st1++) {
                            if (st1 == st2_crf) continue; // not allowed since no change in base
                            if (st2_crf >= NBASE && st1 != (uint8_t)(st2_crf-NBASE)) continue; // can only enter flop base from
                            // same flip base
                            float score = prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1];
                            if (score > curr_score[st2]) {
                                curr_score[st2] = score;
                                traceback[t][st2] = st1;
                            }
                        }
                    }

                    // now look at transitions from end of other words if pos != 0 and pos_in_word = 0
                    if (st2_pos != 0 && st2_pos_in_word == 0) {
                        uint32_t st1_pos = st2_pos - 1;
                        for (uint32_t st1_word_idx = 0; st1_word_idx < wordlen_vocab.size(); st1_word_idx++) {
                            uint32_t st1_pos_in_word = wordlen_vocab[st1_word_idx] - 1;
                            for (uint32_t st1_flip_flop_bit = 0; st1_flip_flop_bit < 2; st1_flip_flop_bit++) {
                                uint32_t st1 = state_tuple2idx[st1_pos][st1_word_idx][st1_pos_in_word][st1_flip_flop_bit];
                                uint8_t st1_crf = vocab[st1_word_idx][st1_pos_in_word] + st1_flip_flop_bit*NBASE;
                                if (st1_crf == st2_crf) continue; // not allowed since no change in base
                                if (st2_crf >= NBASE && st1_crf != st2_crf-NBASE) continue; // can only enter flop base from
                                // same flip base
                                float score = prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
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
                        for (uint32_t st1_flip_flop_bit = 0; st1_flip_flop_bit < 2; st1_flip_flop_bit++) {
                            uint32_t st1 = state_tuple2idx[st1_pos][st1_word_idx][st1_pos_in_word][st1_flip_flop_bit];
                            uint8_t st1_crf = vocab[st1_word_idx][st1_pos_in_word] + st1_flip_flop_bit*NBASE;
                            if (st1_crf == st2_crf) continue; // not allowed since no change in base
                            if (st2_crf >= NBASE && st1_crf != st2_crf-NBASE) continue; // can only enter flop base from
                            // same flip base
                            float score = prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
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
  uint32_t st_pos = nstate_pos-1; // pos is at end for last state
  for (uint32_t st_word_idx = 0; st_word_idx < wordlen_vocab.size(); st_word_idx++) {
      uint32_t st_pos_in_word = wordlen_vocab[st_word_idx] - 1; // last state at end of word
      for (uint32_t st_flip_flop_bit = 0; st_flip_flop_bit < 2; st_flip_flop_bit++) {
          uint32_t st = state_tuple2idx[st_pos][st_word_idx][st_pos_in_word][st_flip_flop_bit];
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
      if (st < nstate_crf) continue; // init state, skip
      auto st_tuple = state_idx2tuple[st];
      if (st_tuple[0] > cur_pos) {
          if (st_tuple[0] != cur_pos + 1) throw std::runtime_error("pos increase not 1");
          if (st_tuple[2] != 0) throw std::runtime_error("pos_in_word at transition not 0");
          cur_pos = st_tuple[0];
          decoded_msg.push_back(st_tuple[1]);
      }
  }
  if(decoded_msg.size() != msg_len) throw std::runtime_error("Decoded message length does not match msg_len");
  return decoded_msg;
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
              float score =
                  prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
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
              if (st2_pos == 0) continue;  // can't go to new state while still being at position 0
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
                  float score =
                      prev_score[st1] + post[t][to_idx_crf_in_post(st2_crf)][st1_crf];
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
