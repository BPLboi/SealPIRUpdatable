#pragma once

#include "pir.hpp"
#include "pir_client.hpp"
#include <map>
#include <memory>
#include <vector>

class PIRServer {
public:
  PIRServer(const seal::EncryptionParameters &enc_params,
            const PirParams &pir_params);

  // NOTE: server takes over ownership of db and frees it when it exits.
  // Caller cannot free db
  void set_database(std::unique_ptr<std::vector<seal::Plaintext>> &&db,
         std::unique_ptr<std::vector<seal::Plaintext>> &&db_non_ntt); // changed from default to also set the non ntt database
  void set_database(const std::unique_ptr<const std::uint8_t[]> &bytes,
                    std::uint64_t ele_num, std::uint64_t ele_size);
  void preprocess_database();

  std::vector<seal::Ciphertext> expand_query(const seal::Ciphertext &encrypted,
                                             std::uint32_t m,
                                             std::uint32_t client_id);

  PirQuery deserialize_query(std::stringstream &stream);
  PirReply generate_reply(PirQuery &query, std::uint32_t client_id);
  // Serializes the reply into the provided stream and returns the number of
  // bytes written
  int serialize_reply(PirReply &reply, std::stringstream &stream);

  void set_galois_key(std::uint32_t client_id, seal::GaloisKeys galkey);

  // Below simple operations are for interacting with the database WITHOUT PIR.
  // So they can be used to modify a particular element in the database or
  // to query a particular element (without privacy guarantees).
  void simple_set(std::uint64_t index, seal::Plaintext pt);
  seal::Ciphertext simple_query(std::uint64_t index);
  void set_one_ct(seal::Ciphertext one);
  
  //----------------------------NEW STUFF:-------------------------------------------------
  seal::Plaintext changed_db_at_idx(const std::unique_ptr<const uint8_t[]> &bytes, std::uint64_t index, std::uint64_t ele_size);
  void update(PirQuery &query, PirReply &reply, std::uint64_t index, seal::Plaintext pt, uint32_t client_id);
  
  
private:
  seal::EncryptionParameters enc_params_; // SEAL parameters
  PirParams pir_params_;                  // PIR parameters
  std::unique_ptr<Database> db_;
  bool is_db_preprocessed_;
  std::map<int, seal::GaloisKeys> galoisKeys_;
  std::unique_ptr<seal::Evaluator> evaluator_;
  std::unique_ptr<seal::BatchEncoder> encoder_;
  std::shared_ptr<seal::SEALContext> context_;

  // This is only used for simple_query
  seal::Ciphertext one_;

  void multiply_power_of_X(const seal::Ciphertext &encrypted,
                           seal::Ciphertext &destination, std::uint32_t index);
  
  //-----------------------NEW STUFF-----------------------------
  std::unique_ptr<Database> non_ntt_db_;
  seal::Ciphertext get_partial_expansion_d1(const seal::Ciphertext &orig,std::uint64_t m, uint64_t index,uint32_t client_id);
};
