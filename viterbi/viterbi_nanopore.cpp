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
#include "cxxopts.hpp"

const uint8_t NBASE = 4;
const char int2base[NBASE] = {'A', 'C', 'G', 'T'};
const uint8_t nstate_crf = 8;
// states are A+,C+,G+,T+,A-,C-,G-,T- (+ = flip, - = flop), you can enter from a
// different base only into a flip base.

typedef std::array<std::array<float, nstate_crf>, NBASE + 1> crf_mat_t;
// for each input state, we have five (4+1) possible output states: e.g., A+ ->
// A+,C+,G+,T+,A-; A- -> A+,C+,G+,T+,A-; A-; etc.

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
// (set in set_conv_params())

uint8_t mem_conv;
uint32_t nstate_conv;
const uint8_t n_out_conv = 2;
uint32_t G[n_out_conv]; // octal
uint32_t initial_state_conv;
uint32_t sync_marker_length;
bool sync_marker[BITSET_SIZE];
uint32_t sync_marker_period;

int set_conv_params(uint8_t mem_conv_param, const std::string sync_marker_param = "", const uint32_t sync_marker_period_param = 0);

std::vector<bool> conv_encode(const std::vector<bool> &msg);

std::vector<bool> read_bit_array(const std::string &infile);

void write_bit_array(const std::vector<bool> &outvec,
                     const std::string &outfile);

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

uint32_t get_state_idx(const uint32_t st_pos, const uint32_t st_conv,
                       const uint32_t st_crf);

uint32_t conv_next_state(const uint32_t cur_state, const bool bit);

uint32_t conv_prev_state(const uint32_t cur_state, const bool bit);

uint32_t conv_output(const uint8_t output_idx, const uint32_t cur_state, const bool bit);

float combine_scores(const float &score_1, const float &score_2, const uint8_t &list_decoding_mode);

std::vector<std::vector<bool>> decode_post_conv_parallel_LVA(
    const std::vector<crf_mat_t> &post,
    const uint32_t msg_len, const uint8_t list_decoding_mode,
    const uint32_t list_size, const uint32_t num_thr);

std::vector<uint32_t> decode_post_vocab(
    const std::vector<crf_mat_t> &post, const uint32_t msg_len,
    const std::vector<std::vector<uint8_t>> &vocab);

int main(int argc, char **argv) {
    cxxopts::Options options("viterbi_nanopore","Viterbi decoder for nanopore dna storage codes");
    options.add_options()
        ("m,mode","Mode: encode_conv, decode_conv, decode_vocab",cxxopts::value<std::string>())
        ("i,infile","Infile with message (encoding) or posterior matrix (decoding)",cxxopts::value<std::string>())
        ("o,outfile","Outfile with encoded/decoded message (list)",cxxopts::value<std::string>())
        ("v,vocabfile","Vocabulary file for decode_vocab mode",cxxopts::value<std::string>())
        ("msg-len","Message length (for decoding)",cxxopts::value<uint32_t>())
        ("mem-conv","Code memory for convolutional code",cxxopts::value<uint8_t>())
        ("sync-marker","Sync marker for convolutional code decoding as string (e.g. 110) (default '')",cxxopts::value<std::string>()->default_value(""))
        ("sync-period","Sync marker period for convolutional code decoding",cxxopts::value<uint32_t>())
        ("l,list-size","List size for convolutional code decoding (default 1)",cxxopts::value<uint32_t>()->default_value("1"))
        ("t,num-thr","Number of threads for convolutional code decoding (default 1)",cxxopts::value<uint32_t>()->default_value("1"))
        ("use-logsumexp","Use log sum exp for combining paths instead of max")
        ("h,help","Display this message")
        ;
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << "\n";
        return 0;
    }
    if (!result.count("mode") ||
            !result.count("infile") ||
            !result.count("outfile")) {
        std::cout << "Invalid options.\n";
        std::cout << options.help() << "\n";
        return -1;
    }
    std::string mode = result["mode"].as<std::string>();
    std::string infile = result["infile"].as<std::string>();
    std::string outfile = result["outfile"].as<std::string>();
    if (mode == "encode_conv" || mode == "decode_conv") {
        if (!result.count("mem-conv")) {
            std::cout << "Memory of convolutional code not specified.\n";
            std::cout << options.help() << "\n";
            return -1;
        }
        int status = 0;
        std::string sync_marker_param = result["sync-marker"].as<std::string>();
        if (sync_marker_param == "")
            status = set_conv_params(result["mem-conv"].as<uint8_t>());
        else
            status = set_conv_params(result["mem-conv"].as<uint8_t>(),sync_marker_param,result["sync-period"].as<uint32_t>());
        if (status != 0) {
            std::cout << options.help() << "\n";
            return -1;
        }
        if (mode == "encode_conv") {
            std::vector<bool> msg = read_bit_array(infile);
            std::vector<bool> encoded_msg = conv_encode(msg);
            write_bit_array_in_bases(encoded_msg, outfile);
        } else {
            // do list decoding
            //
            // modes:
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
               
            if (!result.count("msg-len")) {
                std::cout << "msg-len not specified.\n";
                std::cout << options.help() << "\n";
                return -1;
            }
            std::vector<crf_mat_t> post = read_crf_post(infile);
            uint8_t list_decoding_mode = result["use-logsumexp"].as<bool>()?2:1;
            uint32_t list_size = result["list-size"].as<uint32_t>();
            uint32_t num_thr = result["num-thr"].as<uint32_t>();
            auto decoded_msg_list = decode_post_conv_parallel_LVA(post,result["msg-len"].as<uint32_t>(),list_decoding_mode,list_size,num_thr);
            std::ofstream fout(outfile);
            for (auto decoded_msg : decoded_msg_list) {
                for (auto decoded_msg_bit : decoded_msg)
                    fout << std::to_string(decoded_msg_bit);
                fout << "\n";
            }
        }
    } else if (mode == "decode_vocab") {
        if (!result.count("vocabfile") || !result.count("msg-len")) {
            std::cout << "Vocab file not specified or msg-len not specified.\n";
            std::cout << options.help() << "\n";
            return -1;
        }
        std::string vocabfile = result["vocabfile"].as<std::string>();
        std::vector<crf_mat_t> post = read_crf_post(infile);
        auto vocab = read_vocab_file(vocabfile);
        auto decoded_msg = decode_post_vocab(post, result["msg-len"].as<uint32_t>(), vocab);
        write_vector(decoded_msg, outfile);
        return 0;
    } else {
        std::cout << "Invalid mode.\n";
        std::cout << options.help() << "\n";
        return -1;
    }

  return 0;
}

int set_conv_params(uint8_t mem_conv_param, const std::string sync_marker_param, const uint32_t sync_marker_period_param) {
    mem_conv = mem_conv_param;
    nstate_conv = 1 << mem_conv;
    switch(mem_conv) {
        case 6: G[0] = 0171;
                G[1] = 0133;
                initial_state_conv = 0b100101;
                break;
        case 8: G[0] = 0515;
                G[1] = 0677;
                initial_state_conv = 0b10010110;
                break;
        case 11:G[0] = 05537;
                G[1] = 06131;
                initial_state_conv = 0b10010110001;
                break;
        case 14:G[0] = 075063;
                G[1] = 056711;
                initial_state_conv = 0b10010110001101;
                break;
        default: std::cout << "Invalid mem_conv (allowed: 6, 8, 11, 14)\n";
                      return -1;
    }
    if (sync_marker_param == "")
        return 0;
    
    sync_marker_length = sync_marker_param.length();
    if (sync_marker_length >= BITSET_SIZE) {
        std::cout << "Too large sync marker\n";
        return -1;
    }
    sync_marker_period = sync_marker_period_param;
    if (sync_marker_period < sync_marker_length) {
        std::cout << "Sync period shorter than sync marker length\n";
        return -1;
    }

    for (uint32_t i = 0; i < sync_marker_length; i++) {
        switch(sync_marker_param[i]) {
            case '0': sync_marker[i] = 0;
                      break;
            case '1': sync_marker[i] = 1;
                      break;
            default: std::cout << "Invalid character in sync marker (only 0/1 allowed)\n";
                          return -1;
        }
    }
    return 0;
}

uint32_t conv_next_state(const uint32_t cur_state, const bool bit) {
    if (bit)
        return ((cur_state | nstate_conv) >> 1);
    else
        return (cur_state >> 1);
}

uint32_t conv_prev_state(const uint32_t cur_state, const bool bit) {
    if (bit)
        return ((cur_state << 1) | 1) & (nstate_conv - 1);
    else
        return (cur_state << 1) & (nstate_conv - 1);
}

uint32_t conv_output(const uint8_t output_idx, const uint32_t cur_state, const bool bit) {
    if (bit)
        return __builtin_parity(((cur_state | nstate_conv) & G[output_idx]));
    else
        return __builtin_parity((cur_state & G[output_idx]));
}

std::vector<bool> conv_encode(const std::vector<bool> &msg) {
  std::vector<bool> encoded_msg;
  uint32_t cur_state = initial_state_conv;
  for (bool msg_bit : msg) {
    encoded_msg.push_back(conv_output(0,cur_state,msg_bit));
    encoded_msg.push_back(conv_output(1,cur_state,msg_bit));
    cur_state = conv_next_state(cur_state,msg_bit);
  }
  // add terminating bits
  for (uint8_t i = 0; i < mem_conv; i++) {
    encoded_msg.push_back(conv_output(0,cur_state,0));
    encoded_msg.push_back(conv_output(1,cur_state,0));
    cur_state = conv_next_state(cur_state,0);
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



uint32_t get_state_idx(const uint32_t st_pos, const uint32_t st_conv,
                       const uint32_t st_crf) {
  return st_pos * nstate_conv * nstate_crf + st_conv * nstate_crf + st_crf;
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

std::vector<std::vector<bool>> decode_post_conv_parallel_LVA(
    const std::vector<crf_mat_t> &post,
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
  uint32_t st2_pos_start = std::max((int64_t)nstate_pos - 2 - (nblk - 1 - t), (int64_t)0); 
  uint32_t st2_pos_end = std::min(t + 2, nstate_pos);  
  // only allow pos which can have non -INF scores or will lead to useful
  // final states initially large pos is not allowed, and at the end small
  // pos not allowed (since those can't lead to correct st2_pos at the end).

  for (uint32_t st2_pos = st2_pos_start; st2_pos < st2_pos_end; st2_pos++) {
#pragma omp parallel
#pragma omp for schedule(dynamic)
      for (uint32_t st2_conv = 0; st2_conv < nstate_conv; st2_conv++) {

        // check if this is a valid state, otherwise continue
        bool valid_state = true;
        for (uint32_t shift = 0; shift < mem_conv; shift++) {
            int64_t pos_in_msg = (int64_t)(st2_pos)-1-(int64_t)shift; // can be from -mem_conv to msg_len+mem_conv-1
            bool bit_at_shift = (st2_conv >> (mem_conv - 1 - shift))&1;
            if (pos_in_msg < 0) {
                // initial state should match
                if (bit_at_shift != ((initial_state_conv>> (mem_conv+pos_in_msg))&1)) {
                    valid_state = false;
                    break;
                }
            } else if (pos_in_msg >= msg_len) {
                // final state
                if (bit_at_shift == 1) {
                    valid_state = false;
                    break;
                }
            } else if (sync_marker_length > 0 && pos_in_msg% sync_marker_period < sync_marker_length) {
                // sync marker should match
                if (bit_at_shift != sync_marker[pos_in_msg % sync_marker_period]) {
                    valid_state = false;
                    break;
                }
            }
        }
        if (!valid_state)
            continue;
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

              if (st1_crf < NBASE)
                  endpos_flip[st1_crf] = startpos_flip[st1_crf] = pos_in_Lpl;
              for (uint8_t conv_bit = 0; conv_bit < 2; conv_bit++) {
                uint32_t st1_conv = conv_prev_state(st2_conv,conv_bit);
                if (2 * conv_output(0,st1_conv,curr_conv_bit) +
                        conv_output(1,st1_conv,curr_conv_bit) ==
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


float combine_scores(const float &score_1, const float &score_2, const uint8_t &list_decoding_mode) {
    if (list_decoding_mode == 1) 
        return std::max(score_1,score_2);
    else
        return logsumexpf(score_1,score_2);
}

