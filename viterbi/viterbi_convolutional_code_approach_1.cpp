#include <omp.h>
#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include "cxxopts.hpp"

const uint8_t NBASE = 4;
const char int2base[NBASE] = {'A', 'C', 'G', 'T'};

typedef std::array<float, NBASE + 1> ctc_mat_t;
// at each time step, we have probabilities for A,C,G,T,blank
// uint8_t base: 0->A, 1->C, 2->G, 3->T


// for parallel LVA
const uint32_t BITSET_SIZE = 256;  // 32 bytes
typedef std::bitset<BITSET_SIZE> bitset_t;

float logsumexpf(float x, float y);

// this is the main stucture to store the top paths for each state
struct LVA_path_t {
  bitset_t msg;
  float score_nonblank; // score for path ending with non-blank base
  float score_blank; // score for path ending with non-blank base
  float score; // logsumexp(score_nonblank,score_blank), used for sorting
  // If score is -INF, score_nonblank and score_blank can be garbage 
  uint8_t last_base; // for computing updated score_nonblank for stay transition
                     // and to compare with new_base in non-stay transitions 
                     // for checking whether nonblank->nonblank makes sense

  LVA_path_t() {
    float INF = std::numeric_limits<float>::infinity();
    score_nonblank = -INF;
    score_blank = -INF;
    compute_score();
  }

  LVA_path_t(const bitset_t &msg_, const float &score_nonblank_, 
             const float &score_blank_, const uint8_t &last_base_) {
    msg = msg_;
    score_nonblank = score_nonblank_;
    score_blank = score_blank_;
    last_base = last_base_;
  }

  // update score based on nonblank and blank score.
  // NOTE: this must be called externally
  void compute_score() {
    score = logsumexpf(score_nonblank, score_blank);
  }
};

// struct for storing information about previous state for a current state and
// the transition
struct prev_state_info_t {
  uint32_t st_conv;
  uint8_t new_base; // undefined for transition from same conv state, otherwise
                    // stores the new base added in this transition 
  uint8_t msg_shift;  // shift in message in transition
  uint8_t msg_newbits;  // new bits in transition (new_msg = (old_msg <<
                        // msg_shift)|msg_newbits)
};

bool rc_flag = false;

// convolutional code related parameters
// (set in set_conv_params())

uint8_t mem_conv;
uint32_t nstate_conv;
const uint8_t n_out_conv = 2;
uint32_t G[n_out_conv];  // octal
uint32_t initial_state_conv;
uint32_t final_state_conv;
uint32_t sync_marker_length;
bool sync_marker[BITSET_SIZE];
uint32_t sync_marker_period;
const uint8_t MAX_PUNCTURING_PATTERN_SIZE = 16;
uint8_t puncturing_pattern[MAX_PUNCTURING_PATTERN_SIZE];
uint8_t puncturing_pattern_len;
uint32_t nstate_pos;
uint32_t st_pos2msg_pos[BITSET_SIZE];  // map st_pos to position in message

int set_conv_params(uint8_t mem_conv_param, uint8_t rate_param,
                    uint32_t msg_len_param,
                    const std::string sync_marker_param = "",
                    const uint32_t sync_marker_period_param = 0);

std::vector<bool> conv_encode(const std::vector<bool> &msg);

std::vector<std::vector<bool>> read_bit_array(const std::string &infile);

void write_bit_array(const std::vector<bool> &outvec,
                     const std::string &outfile);

void write_bit_array_in_bases(const std::vector<std::vector<bool>> &outvec_vec,
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
  float INF = std::numeric_limits<float>::infinity();
  if (x == -INF && y == -INF)
    return -INF; // the formula below returns nan
  float max_x_y = std::max(x, y);
  return max_x_y + logf(expf(x - max_x_y) + expf(y - max_x_y));
}

std::vector<ctc_mat_t> read_ctc_post(const std::string &infile);

std::vector<std::vector<uint8_t>> read_vocab_file(const std::string &infile);

uint32_t conv_next_state(const uint32_t cur_state, const bool bit);

uint32_t conv_prev_state(const uint32_t cur_state, const bool bit);

uint32_t conv_output(const uint8_t output_idx, const uint32_t cur_state,
                     const bool bit);

bool is_valid_state(const uint32_t &st2_pos, const uint32_t &st2_conv,
                    const uint32_t &msg_len);

std::vector<std::vector<bool>> decode_post_conv_parallel_LVA(
    const std::vector<ctc_mat_t> &post, const uint32_t msg_len,
    const uint32_t list_size, const uint32_t num_thr,
    const uint32_t max_deviation);

std::vector<prev_state_info_t> find_prev_states(const uint32_t &st2_conv,
                                                const uint8_t &punc_pattern);

uint32_t reverse_integer_bits(const uint32_t &num, const uint32_t numbits);

int main(int argc, char **argv) {
  cxxopts::Options options("viterbi_nanopore",
                           "Viterbi decoder for nanopore dna storage codes");
  options.add_options()("m,mode", "Mode: encode, decode",
                        cxxopts::value<std::string>())(
      "i,infile",
      "Infile with message (encoding) or posterior matrix (decoding)",
      cxxopts::value<std::string>())(
      "o,outfile", "Outfile with encoded/decoded message (list)",
      cxxopts::value<std::string>())("msg-len", "Message length",
                                     cxxopts::value<uint32_t>())(
      "mem-conv", "Code memory for convolutional code",
      cxxopts::value<uint8_t>())(
      "sync-marker",
      "Sync marker for convolutional code decoding as string (e.g. 110) "
      "(default '')",
      cxxopts::value<std::string>()->default_value(""))(
      "sync-period", "Sync marker period for convolutional code decoding",
      cxxopts::value<uint32_t>())(
      "l,list-size", "List size for convolutional code decoding (default 1)",
      cxxopts::value<uint32_t>()->default_value("1"))(
      "r,rate",
      "Rate of convolutional code: options 1 (1/2), 2 (2/3), 3 (3/4), 4 (4/5), "
      "5 (5/6), 7 (7/8) "
      "(default 1). Use standard puncturing patterns, expects appropriate "
      "padding (at most 1 bit needed) to make output length even.",
      cxxopts::value<uint8_t>()->default_value("1"))(
      "max-deviation",
      "Max allowable deviation of st_pos around its expected value during "
      "decoding (tradeoff b/w speed and accuracy) (default: infinite)",
      cxxopts::value<uint32_t>())("rc",
                                  "Reverse complement read (for decoding)")(
      "t,num-thr",
      "Number of threads for convolutional code decoding (default 1)",
      cxxopts::value<uint32_t>()->default_value("1"))("h,help",
                                                      "Display this message");
  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << "\n";
    return 0;
  }
  if (!result.count("mode") || !result.count("infile") ||
      !result.count("outfile")) {
    std::cout << "Invalid options.\n";
    std::cout << options.help() << "\n";
    return -1;
  }
  std::string mode = result["mode"].as<std::string>();
  std::string infile = result["infile"].as<std::string>();
  std::string outfile = result["outfile"].as<std::string>();
  if (mode == "encode" || mode == "decode") {
    if (!result.count("mem-conv")) {
      std::cout << "Memory of convolutional code not specified.\n";
      std::cout << options.help() << "\n";
      return -1;
    }
    if (!result.count("msg-len")) {
      std::cout << "msg-len not specified.\n";
      std::cout << options.help() << "\n";
      return -1;
    }
    uint32_t msg_len = result["msg-len"].as<uint32_t>();
    int status = 0;
    if (mode == "decode") rc_flag = result["rc"].as<bool>();
    if (rc_flag) std::cout << "Reverse complement flag detected.\n";
    std::string sync_marker_param = result["sync-marker"].as<std::string>();
    if (sync_marker_param == "")
      status = set_conv_params(result["mem-conv"].as<uint8_t>(),
                               result["rate"].as<uint8_t>(), msg_len);
    else
      status = set_conv_params(
          result["mem-conv"].as<uint8_t>(), result["rate"].as<uint8_t>(),
          msg_len, sync_marker_param, result["sync-period"].as<uint32_t>());
    if (status != 0) {
      std::cout << options.help() << "\n";
      return -1;
    }
    if (mode == "encode") {
      std::vector<std::vector<bool>> msg_vec = read_bit_array(infile);
      std::vector<std::vector<bool>> encoded_msg_vec;
      for (auto msg : msg_vec) {
        if (msg.size() != msg_len) {
          std::cout << "Message length does not match msg_len parameter.\n";
          return -1;
        }
        encoded_msg_vec.push_back(conv_encode(msg));
      }
      write_bit_array_in_bases(encoded_msg_vec, outfile);
    } else {
      // do list decoding
      // 
      // based on ideas from Parallel LVA algorithm as described in 
      // https://github.com/shubhamchandak94/kBestViterbi/blob/master/kBestViterbi.py
      // or in ieeexplore.ieee.org/iel1/26/12514/00577040.pdf
      // and on beam search for CTC as described in https://distill.pub/2017/ctc/
      // and in https://gist.github.com/awni/56369a90d03953e370f3964c826ed4b0
      // 
      // We have states correspoding to conv code state and pos in msg.
      // Each state stores a list of messages with their score of ending in 
      // blank and non-blank. At next step, we add one character, add (logsumexp) the 
      // scores for all the ways a new message can be obtained, and take the top ones.

      uint32_t max_deviation =
          msg_len + mem_conv +
          1;  // don't restrict anything, do full exact Viterbi
      if (result.count("max-deviation"))
        max_deviation = result["max-deviation"].as<uint32_t>();
      std::vector<ctc_mat_t> post = read_ctc_post(infile);
      uint32_t list_size = result["list-size"].as<uint32_t>();
      uint32_t num_thr = result["num-thr"].as<uint32_t>();
      auto decoded_msg_list = decode_post_conv_parallel_LVA(
          post, msg_len, list_size, num_thr, max_deviation);
      std::ofstream fout(outfile);
      for (auto decoded_msg : decoded_msg_list) {
        for (auto decoded_msg_bit : decoded_msg)
          fout << std::to_string(decoded_msg_bit);
        fout << "\n";
      }
    }
  } else {
    std::cout << "Invalid mode.\n";
    std::cout << options.help() << "\n";
    return -1;
  }

  return 0;
}

int set_conv_params(uint8_t mem_conv_param, uint8_t rate_param,
                    uint32_t msg_len_param, const std::string sync_marker_param,
                    const uint32_t sync_marker_period_param) {
  mem_conv = mem_conv_param;
  nstate_conv = 1 << mem_conv;
  switch (mem_conv) {
    case 6:
      G[0] = 0171;
      G[1] = 0133;
      initial_state_conv = 0b100101;
      break;
    case 8:
      G[0] = 0515;
      G[1] = 0677;
      initial_state_conv = 0b10010110;
      break;
    case 11:
      G[0] = 05537;
      G[1] = 06131;
      initial_state_conv = 0b10010110001;
      break;
    case 14:
      G[0] = 075063;
      G[1] = 056711;
      initial_state_conv = 0b10010110001101;
      break;
    default:
      std::cout << "Invalid mem_conv (allowed: 6, 8, 11, 14)\n";
      return -1;
  }
  final_state_conv = reverse_integer_bits(initial_state_conv, mem_conv);

  // puncturing pattern representing using an array: 0 representing 1\n1, 1
  // representing 0,1\n1,0, 2 representing 1,0\n0,1, 3 representing 0,0\n1,1
  // (these are the only patterns supported right now)
  switch (rate_param) {
    case 1:
      puncturing_pattern[0] = 0;
      puncturing_pattern_len = 1;
      break;
    case 2:
      puncturing_pattern[0] = 0;
      puncturing_pattern[1] = 2;
      puncturing_pattern[2] = 0;
      puncturing_pattern_len = 3;
      break;
    case 3:
      puncturing_pattern[0] = 0;
      puncturing_pattern[1] = 1;
      puncturing_pattern_len = 2;
      break;
    case 4:
      puncturing_pattern[0] = 0;
      puncturing_pattern[1] = 3;
      puncturing_pattern[2] = 0;
      puncturing_pattern[3] = 2;
      puncturing_pattern[4] = 1;
      puncturing_pattern_len = 5;
      break;
    case 5:
      puncturing_pattern[0] = 0;
      puncturing_pattern[1] = 1;
      puncturing_pattern[2] = 2;
      puncturing_pattern_len = 3;
      break;
    case 7:
      puncturing_pattern[0] = 0;
      puncturing_pattern[1] = 3;
      puncturing_pattern[2] = 1;
      puncturing_pattern[3] = 1;
      puncturing_pattern_len = 4;
      break;
    default:
      std::cout << "Invalid rate parameter (allowed: 1, 3, 5, 7)\n";
      return -1;
  }
  // verify that msg_len is such that output length is even (for converting to
  // DNA). also set nstate_pos (which corresponds now to the output position
  // rather than input). also set st_pos2msg_pos (where the msg_pos is indexed
  // by 1 (so ranges from 1 to msg_len_param + mem_conv)).
  nstate_pos = 1;  // state 0 (init)
  st_pos2msg_pos[0] = 0;
  uint32_t i = 0;
  while (i < msg_len_param + mem_conv) {
    i += ((puncturing_pattern[(nstate_pos - 1) % puncturing_pattern_len] == 0)
              ? 1
              : 2);
    st_pos2msg_pos[nstate_pos++] = i;
  }
  if (i != msg_len_param + mem_conv) {
    std::cout
        << "Output length not even. Try padding with a single 0 at end.\n";
    return -1;
  }

  if (rc_flag) {
    // basically reverse everything
    G[0] = reverse_integer_bits(G[0], mem_conv + 1);
    G[1] = reverse_integer_bits(G[1], mem_conv + 1);
    uint32_t temp = reverse_integer_bits(initial_state_conv, mem_conv);
    initial_state_conv = reverse_integer_bits(final_state_conv, mem_conv);
    final_state_conv = temp;

    // now do st_pos2msg_pos and puncturing_pattern
    uint32_t temp_puncturing_pattern[BITSET_SIZE];
    for (uint32_t i = 0; i < puncturing_pattern_len; i++)
      temp_puncturing_pattern[i] = puncturing_pattern[i];
    // first find where end of forward read lies in terms of puncturing pattern
    uint32_t punc_pattern_end_idx =
        (nstate_pos - 2) %
        puncturing_pattern_len;  // -2 because st_pos is 1 indexed
    uint8_t punc_pattern_reverse_map[4] = {
        0, 2, 1, 3};  // 1 <-> 2 when pattern looked in opposite direction
    for (uint32_t i = 0; i < puncturing_pattern_len; i++) {
      uint32_t j = (puncturing_pattern_len - i + punc_pattern_end_idx) %
                   puncturing_pattern_len;
      puncturing_pattern[i] =
          punc_pattern_reverse_map[temp_puncturing_pattern[j]];
    }
    std::reverse(st_pos2msg_pos, st_pos2msg_pos + nstate_pos);
    for (uint32_t i = 0; i < nstate_pos; i++)
      st_pos2msg_pos[i] = msg_len_param + mem_conv - st_pos2msg_pos[i];
  }

  if (sync_marker_param == "") return 0;

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
    switch (sync_marker_param[i]) {
      case '0':
        sync_marker[i] = 0;
        break;
      case '1':
        sync_marker[i] = 1;
        break;
      default:
        std::cout << "Invalid character in sync marker (only 0/1 allowed)\n";
        return -1;
    }
  }
  return 0;
}

uint32_t reverse_integer_bits(const uint32_t &num, const uint32_t numbits) {
  uint32_t ret = 0;
  for (uint32_t i = 0; i < numbits; i++) {
    ret <<= 1;
    ret |= ((num >> i) & 1);
  }
  return ret;
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

uint32_t conv_output(const uint8_t output_idx, const uint32_t cur_state,
                     const bool bit) {
  // XOR with rc_flag to capture reverse complementation of base
  if (bit)
    return (__builtin_parity(((cur_state | nstate_conv) & G[output_idx]))) ^
           (rc_flag);
  else
    return (__builtin_parity((cur_state & G[output_idx]))) ^ (rc_flag);
}

std::vector<bool> conv_encode(const std::vector<bool> &msg) {
  std::vector<bool> encoded_msg;
  uint32_t cur_state = initial_state_conv;
  for (bool msg_bit : msg) {
    encoded_msg.push_back(conv_output(0, cur_state, msg_bit));
    encoded_msg.push_back(conv_output(1, cur_state, msg_bit));
    cur_state = conv_next_state(cur_state, msg_bit);
  }
  // add terminating bits
  for (uint8_t i = 0; i < mem_conv; i++) {
    bool padding_bit = (final_state_conv >> i) & 1;
    encoded_msg.push_back(conv_output(0, cur_state, padding_bit));
    encoded_msg.push_back(conv_output(1, cur_state, padding_bit));
    cur_state = conv_next_state(cur_state, padding_bit);
  }
  if (cur_state != final_state_conv)
    throw std::runtime_error(
        "state after encoding does not match final_state_conv");
  std::vector<bool> encoded_msg_punctured;
  // puncture
  uint32_t i = 0;  // position in encoded message
  for (uint32_t st_pos = 0; st_pos < nstate_pos - 1; st_pos++) {
    uint8_t punc_pattern = puncturing_pattern[st_pos % puncturing_pattern_len];
    switch (punc_pattern) {
      case 0:
        encoded_msg_punctured.push_back(encoded_msg[i]);
        encoded_msg_punctured.push_back(encoded_msg[i + 1]);
        i += 2;
        break;
      case 1:
        encoded_msg_punctured.push_back(encoded_msg[i + 1]);
        encoded_msg_punctured.push_back(encoded_msg[i + 2]);
        i += 4;
        break;
      case 2:
        encoded_msg_punctured.push_back(encoded_msg[i]);
        encoded_msg_punctured.push_back(encoded_msg[i + 3]);
        i += 4;
        break;
      case 3:
        encoded_msg_punctured.push_back(encoded_msg[i + 1]);
        encoded_msg_punctured.push_back(encoded_msg[i + 3]);
        i += 4;
        break;
    }
  }
  if (i != encoded_msg.size())
    throw std::runtime_error("error in encoding, lengths don't match");
  return encoded_msg_punctured;
}

std::vector<std::vector<bool>> read_bit_array(const std::string &infile) {
  std::ifstream fin(infile);
  std::vector<std::vector<bool>> vec_vec;
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
      case '\n':
        vec_vec.push_back(vec);
        vec.clear();
        break;
      default:
        throw std::runtime_error("invalid character in input file");
    }
  }
  fin.close();
  return vec_vec;
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

void write_bit_array_in_bases(const std::vector<std::vector<bool>> &outvec_vec,
                              const std::string &outfile) {
  std::ofstream fout(outfile);
  for (auto outvec : outvec_vec) {
    uint32_t len = outvec.size();
    if (len % 2 != 0) throw std::runtime_error("length not even");
    for (uint32_t i = 0; i < len / 2; i++)
      fout << int2base[2 * outvec[2 * i] + outvec[2 * i + 1]];
    fout << "\n";
  }
  fout.close();
}

std::vector<ctc_mat_t> read_ctc_post(const std::string &infile) {
  // Read post file - note that post file has probabilities
  // so we convert to log probabilities
  std::ifstream fin(infile, std::ios::binary);
  std::vector<ctc_mat_t> post;
  ctc_mat_t post_mat;
  float val;
  fin.read((char *)&val, sizeof(float));
  while (!fin.eof()) {
    post_mat[NBASE] = logf(val);
    fin.read((char *)&val, sizeof(float));
    for (uint8_t i = 0; i < NBASE; i++) {
      post_mat[i] = logf(val);
      fin.read((char *)&val, sizeof(float));
    }
    post.push_back(post_mat);
  }
  return post;
}

uint32_t get_state_idx(const uint32_t st_pos, const uint32_t st_conv) {
  return st_pos * nstate_conv + st_conv;
}

std::vector<std::vector<bool>> decode_post_conv_parallel_LVA(
    const std::vector<ctc_mat_t> &post, const uint32_t msg_len,
    const uint32_t list_size, const uint32_t num_thr,
    const uint32_t max_deviation) {
  omp_set_num_threads(num_thr);
  float INF = std::numeric_limits<float>::infinity();
  uint64_t nstate_total_64 = nstate_pos * nstate_conv;
  if (nstate_total_64 >= ((uint64_t)1 << 32))
    throw std::runtime_error("Too many states, can't fit in 32 bits");
  uint32_t nstate_total = (uint32_t)nstate_total_64;
  uint32_t nblk = post.size();
  if (post.size() < nstate_pos + 1)
    throw std::runtime_error("Too small post matrix");

  // instead of traceback, store the msg till now as a bitset
  if (msg_len > BITSET_SIZE)
    throw std::runtime_error("msg_len can't be above BITSET_SIZE");

  // arrays for storing previous and current best paths
  // [nstate_total][list_size]
  LVA_path_t **curr_best_paths = new LVA_path_t *[nstate_total];
  LVA_path_t **prev_best_paths = new LVA_path_t *[nstate_total];
  for (uint32_t i = 0; i < nstate_total; i++) {
    curr_best_paths[i] = new LVA_path_t[list_size]();
    prev_best_paths[i] = new LVA_path_t[list_size]();
  }

  // lambda expression to compare paths (decreasing in score)
  auto LVA_path_t_compare = [](const LVA_path_t &a, const LVA_path_t &b) {
                              return a.score > b.score;
                            };

  // find valid states based on intial and final states as well as
  // synchronization markers
  std::vector<bool> valid_state_array(nstate_pos * nstate_conv);
#pragma omp parallel
#pragma omp for
  for (uint32_t st_pos = 0; st_pos < nstate_pos; st_pos++)
    for (uint32_t st_conv = 0; st_conv < nstate_conv; st_conv++)
      valid_state_array[nstate_conv * st_pos + st_conv] =
          is_valid_state(st_pos2msg_pos[st_pos], st_conv, msg_len);

  // precompute the previous states and associated info for all states now
  // note that this is valid only for st_pos > 0 (if st_pos = 0, only previous
  // state allowed is same state - which is always first entry in the
  // prev_state_vector)
  std::vector<std::vector<std::vector<prev_state_info_t>>> prev_state_vector(4);
#pragma omp parallel
#pragma omp for
  for (uint8_t punc_pattern = 0; punc_pattern < 4; punc_pattern++) {
    // only fill if this punc_pattern is relevant for this rate
    if (std::find(puncturing_pattern,
                  puncturing_pattern + puncturing_pattern_len,
                  punc_pattern) == puncturing_pattern + puncturing_pattern_len)
      continue;
    prev_state_vector[punc_pattern].resize(nstate_conv);
    for (uint32_t st_conv = 0; st_conv < nstate_conv; st_conv++)
        prev_state_vector[punc_pattern][st_conv] =
            find_prev_states(st_conv, punc_pattern);
  }

  // set score_blank to zero for initial state
  uint32_t initial_st = get_state_idx(0, initial_state_conv);
  curr_best_paths[initial_st][0].score_blank = 0.0;
  curr_best_paths[initial_st][0].score_nonblank = -INF;
  curr_best_paths[initial_st][0].compute_score();
  // forward Viterbi pass

  for (uint32_t t = 0; t < nblk; t++) {
    // swap prev and curr arrays
    std::swap(curr_best_paths, prev_best_paths);

    // st is current state
    uint32_t st_pos_start =
        std::max((int64_t)nstate_pos - 2 - (nblk - 1 - t), (int64_t)0);
    uint32_t st_pos_end = std::min(t + 2, nstate_pos);

    st_pos_start = std::max(
        (int64_t)st_pos_start, (int64_t)((double)(t) / nblk * nstate_pos - max_deviation));
    st_pos_end = std::min(st_pos_start + 2 * max_deviation, st_pos_end);

    // only allow pos which can have non -INF scores or will lead to useful
    // final states initially large pos is not allowed, and at the end small
    // pos not allowed (since those can't lead to correct st_pos at the end).

#pragma omp parallel
#pragma omp for schedule(dynamic)
    for (uint32_t st_pos = st_pos_start; st_pos < st_pos_end; st_pos++) {

      // vector containing the candidate items for next step list
      std::vector<LVA_path_t> candidate_paths;

      uint8_t punc_pattern = 0;
      if (st_pos != 0)
        punc_pattern =
            puncturing_pattern[(st_pos - 1) % puncturing_pattern_len];

      for (uint32_t st_conv = 0; st_conv < nstate_conv; st_conv++) {
        // check if this is a valid state, otherwise continue
        if (!valid_state_array[nstate_conv * st_pos + st_conv]) continue;

        // clear candidate_paths
        candidate_paths.clear();

        // first do stay transition
        uint32_t st = get_state_idx(st_pos, st_conv);
        uint32_t num_stay_candidates = 0;
        for (uint32_t i = 0; i < list_size; i++) {
          if (prev_best_paths[st][i].score == -INF)
            break;
          num_stay_candidates++;
          float new_score_blank = logsumexpf(prev_best_paths[st][i].score_blank + post[t][NBASE], 
                                   prev_best_paths[st][i].score_nonblank + post[t][NBASE]);
          float new_score_nonblank = prev_best_paths[st][i].score_nonblank +
                                     post[t][prev_best_paths[st][i].last_base];
          candidate_paths.emplace_back(prev_best_paths[st][i].msg,new_score_nonblank,
                                       new_score_blank,prev_best_paths[st][i].last_base);
        }

        // Now go through non-stay transitions.
        // For each case, first look in the stay transitions to check if the msg has
        // already appeared before
        // start with psidx = 1, since 0 corresponds to stay (already done above)
        const auto &prev_states_st =
            prev_state_vector[punc_pattern][st_conv];
        if (st_pos != 0) { // otherwise only stay transition makes sense
          for (uint32_t psidx = 1; psidx < prev_states_st.size(); psidx++) {
            uint32_t prev_st_pos = st_pos - 1;
            uint32_t prev_st = get_state_idx(prev_st_pos, prev_states_st[psidx].st_conv);
            uint8_t msg_shift = prev_states_st[psidx].msg_shift;
            uint8_t msg_newbits = prev_states_st[psidx].msg_newbits;
            uint8_t new_base = prev_states_st[psidx].new_base;
            for (uint32_t i = 0; i < list_size; i++) {
              if (prev_best_paths[prev_st][i].score == -INF)
                break;
              bitset_t msg = prev_best_paths[prev_st][i].msg;
              msg = (msg << msg_shift) | bitset_t(msg_newbits);
              float new_score_blank = -INF;
              float new_score_nonblank;
              if (new_base != prev_best_paths[prev_st][i].last_base) {
                new_score_nonblank = 
                    logsumexpf(prev_best_paths[prev_st][i].score_blank + post[t][new_base],
                             prev_best_paths[prev_st][i].score_nonblank + post[t][new_base]);
              } else {
                // the newly added base is same as last base so we can't have the 
                // thing ending with non_blank (otherwise it gets collapsed)
                new_score_nonblank = prev_best_paths[prev_st][i].score_blank + post[t][new_base];
              }
              if (new_score_nonblank == -INF)
              {
                continue;
                // overall score is -INF (this might happen if 
                // score_blank for previous path is -INF and we are
                // in second case above) 
              }

              // now check if this is already present in the stay transitions.
              // first try to match the new_base for speed, then look at full msg.
              // if already present, update the nonblank score
              bool match_found = false;
              for (uint32_t j = 0; j < num_stay_candidates; j++) {
                if (new_base == candidate_paths[j].last_base) {
                  if (msg == candidate_paths[j].msg) {
                    match_found = true;
                    candidate_paths[j].score_nonblank = 
                            logsumexpf(candidate_paths[j].score_nonblank,new_score_nonblank);
                    break;
                  }
                }
              }
              if (!match_found) {
                candidate_paths.emplace_back(msg,new_score_nonblank,new_score_blank,new_base);
              }
            }
          }
        }
       
        uint32_t num_candidates = candidate_paths.size(); 
        // update scores based on score_blank and score_nonblank
        for (uint32_t i = 0; i < num_candidates; i++)
           candidate_paths[i].compute_score();

        auto num_candidates_to_keep = std::min(list_size, num_candidates);
        // use nth_element to to do partial sorting if num_candidates_to_keep < num_candidates
        if (num_candidates_to_keep < num_candidates && num_candidates_to_keep > 0)
          std::nth_element(candidate_paths.begin(),
                           candidate_paths.begin()+num_candidates_to_keep-1,
                           candidate_paths.end(),
                           LVA_path_t_compare);
        // copy over top candidate paths to curr_best
        std::copy(candidate_paths.begin(),candidate_paths.begin()+num_candidates_to_keep,
                    curr_best_paths[st]);
        // fill any remaining positions in list with score -INF so they are not used later
        for (uint32_t i = num_candidates_to_keep; i < list_size; i++)
          curr_best_paths[st][i].score = -INF;
      }
    }
  }

  uint32_t st_pos = nstate_pos - 1, st_conv = final_state_conv;  // last state
  LVA_path_t *LVA_path_list_final = curr_best_paths[get_state_idx(st_pos, st_conv)];

  // sort in decreasing order by score 
  // NOTE: the curr_best_paths list is not sorted since we use nth_element partial sorting
  std::sort(LVA_path_list_final, LVA_path_list_final+list_size, LVA_path_t_compare);

  std::vector<std::vector<bool>> decoded_msg_list;

  // now convert bitset to bool vectors
  for (uint32_t list_pos = 0; list_pos < list_size; list_pos++) {
    std::vector<bool> decoded_msg(msg_len);
    for (uint8_t i = 0; i < msg_len; i++)
      decoded_msg[i] =
          LVA_path_list_final[list_pos].msg[msg_len + mem_conv - 1 -
                       i];  // due to way bitset is stored in reverse
    if (rc_flag) std::reverse(decoded_msg.begin(), decoded_msg.end());
    decoded_msg_list.push_back(decoded_msg);
    // FOR DEBUGGING
    /*
        std::cout << "score: " << LVA_path_list_final[list_pos].score << "\n";
        for (auto b : decoded_msg_list.back()) std::cout << b;
        std::cout << "\n\n";
    */
  }
  // std::cout << "Final list size: " << decoded_msg_list.size() << "\n";

  for (uint32_t i = 0; i < nstate_total; i++) {
    delete[] curr_best_paths[i];
    delete[] prev_best_paths[i];
  }
  delete[] curr_best_paths;
  delete[] prev_best_paths;
  return decoded_msg_list;
}

std::vector<prev_state_info_t> find_prev_states(const uint32_t &st2_conv,
                                                const uint8_t &punc_pattern) {
  std::vector<prev_state_info_t> prev_state_vec;
  prev_state_info_t prev_state_info;
  uint8_t curr_conv_bit = (st2_conv >> (mem_conv - 1));
  uint8_t curr_conv_bit_1 = (st2_conv >> (mem_conv - 2)) & 1;
  // first do stay
  uint32_t st1_conv = st2_conv;
  prev_state_info.st_conv = st1_conv;
  prev_state_info.msg_shift = 0;
  prev_state_info.msg_newbits = 0;
  prev_state_vec.push_back(prev_state_info);

  // now depending on the puncturing pattern, find the previous conv_state
  // and the base outputted on the corresponding transition
  if (punc_pattern == 0) {
    for (uint8_t conv_bit = 0; conv_bit < 2; conv_bit++) {
      uint32_t st1_conv = conv_prev_state(st2_conv, conv_bit);
      uint8_t new_base = 2 * conv_output(0, st1_conv, curr_conv_bit)
                         + conv_output(1, st1_conv, curr_conv_bit);
      prev_state_info.st_conv = st1_conv;
      prev_state_info.new_base = new_base;
      prev_state_info.msg_shift = 1;
      prev_state_info.msg_newbits = curr_conv_bit;
      prev_state_vec.push_back(prev_state_info);
    }
  } else {
    for (uint8_t conv_bit = 0; conv_bit < 2; conv_bit++) {
      for (uint8_t conv_bit_1 = 0; conv_bit_1 < 2; conv_bit_1++) {
        uint32_t st1_5_conv = conv_prev_state(st2_conv, conv_bit);
        uint32_t st1_conv = conv_prev_state(st1_5_conv, conv_bit_1);
        uint8_t bit_0, bit_1, bit_2, bit_3;
        bit_0 = conv_output(0, st1_conv, curr_conv_bit_1);
        bit_1 = conv_output(1, st1_conv, curr_conv_bit_1);
        bit_2 = conv_output(0, st1_5_conv, curr_conv_bit);
        bit_3 = conv_output(1, st1_5_conv, curr_conv_bit);
        uint8_t base = 0;
        switch (punc_pattern) {
          case 1:
            base = rc_flag ? (2 * bit_2 + bit_1) : (2 * bit_1 + bit_2);
            break;
          case 2:
            base = rc_flag ? (2 * bit_3 + bit_0) : (2 * bit_0 + bit_3);
            break;
          case 3:
            base = rc_flag ? (2 * bit_3 + bit_1) : (2 * bit_1 + bit_3);
            break;
        }
        prev_state_info.st_conv = st1_conv;
        prev_state_info.new_base = base;
        prev_state_info.msg_shift = 2;
        prev_state_info.msg_newbits = 2 * curr_conv_bit_1 + curr_conv_bit;
        prev_state_vec.push_back(prev_state_info);
      }
    }
  }
  prev_state_vec.shrink_to_fit();
  return prev_state_vec;
}

bool is_valid_state(const uint32_t &msg_pos, const uint32_t &st_conv,
                    const uint32_t &msg_len) {
  bool valid_state = true;
  for (uint32_t shift = 0; shift < mem_conv; shift++) {
    int64_t pos_in_msg =
        (int64_t)(msg_pos)-1 -
        (int64_t)shift;  // can be from -mem_conv to msg_len+mem_conv-1
    int64_t pos_in_msg_fwd = pos_in_msg;  // for sync marker, convert to forward
                                          // strand position if rc_flag
    if (rc_flag) pos_in_msg_fwd = msg_len - 1 - pos_in_msg;
    bool bit_at_shift = (st_conv >> (mem_conv - 1 - shift)) & 1;
    if (pos_in_msg < 0) {
      // initial state should match
      if (bit_at_shift !=
          ((initial_state_conv >> (mem_conv + pos_in_msg)) & 1)) {
        valid_state = false;
        break;
      }
    } else if (pos_in_msg >= msg_len) {
      // final state
      if (bit_at_shift != ((final_state_conv >> (pos_in_msg - msg_len)) & 1)) {
        valid_state = false;
        break;
      }
    } else if (sync_marker_length > 0 &&
               pos_in_msg_fwd % sync_marker_period < sync_marker_length) {
      // sync marker should match
      if (bit_at_shift != sync_marker[pos_in_msg_fwd % sync_marker_period]) {
        valid_state = false;
        break;
      }
    }
  }
  return valid_state;
}
