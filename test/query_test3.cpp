#include "pir.hpp"
#include "pir_client.hpp"
#include "pir_server.hpp"
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <cmath>
#include <seal/seal.h>

using namespace std::chrono;
using namespace std;
using namespace seal;

int query_test_w_update(uint64_t num_items, uint64_t item_size, uint32_t degree,
               uint32_t lt, uint32_t dim);
int getNonCheck(); // gets a random non-check index
int getCheck();
uint64_t sumAllTimes (int trialNum);
double allResultsSTDEV(int trialNum);
void printWhiskerPlot(int trialNum);
// --------------------------------VARIABLES TO CHANGE------------------------------------
static const int databaseSize = 1<< 16; 
static const int checks = 5; // number of database indexes for which the answer will be re-computed each time and re-checked
static const int numUpdates = 4; //number of times a database entry will be updated. half of these will be on random database entries that should cause no chantge, half will be on one of the database indexes that are being changed
static const int numTrials = 1; //number of times all of this will be re-run

//--------------------------------------------------------------------------------------
int checkIdxs[checks];
int allResults[numTrials][numUpdates*checks];
int trialNum = 0;

int main(int argc, char *argv[]) {
  //(query_test(1 << 20, 288, 4096, 20, 2) == 0); 
  /*query_test takes in number of mailboxes, number of bytes that mailboxes will contain, degree of BFV polynomials, log of plaintext modulus, recursion level
  Haven't played around with changing degree of BFV polynomials of the log of the plaintext modulus - just using default values
  For recursion level, the github page says that if the number of BFV plaintexts (= degree of BFV polynomials * log of plaintext modulus / bytes per mailbox) is greater than the degree of the BFV polynomials, d=2 is optimal. Otherwise, d=1 minimizes communication costs, though d=2 gives faster computation costs (faster expand time).
  
  Currently using d=2 due to potentially faster computation costs, though the database sizes mean d=1 would give faster communication costs. Need to do more testing to see if using d=2 is actually a good choice.
  */
  //generates check indexes
  for(int i = 0; i< checks; i++){
    int n;
    while( (n = rand()) > RAND_MAX - RAND_MAX%databaseSize){}
    n = n%databaseSize;
    checkIdxs[i] = n;
    for(int j = 0; j<i; j++){
      if(checkIdxs[j] == checkIdxs[i]) i--;
    }
  }
  
  for(int i = 0; i< checks; i++){
    cout << checkIdxs[i] << endl;
  }
  
  for(int i = 0; i< numTrials; i++){
    assert(query_test_w_update(databaseSize, 288, 4096, 20, 1) == 0);
  }
}

void printWhiskerPlot(int trialNum){
  int n = numUpdates*checks; // stores number of elements in the array
  	// to make things easier to type
  sort(allResults[trialNum],allResults[trialNum] + n);
  
  uint32_t min = allResults[trialNum][0];
  uint32_t max = allResults[trialNum][n-1];
  
  uint32_t median;
  uint32_t q1; 
  uint32_t q3;
  
  if(n%2 == 0){
   median = (allResults[trialNum][n/2] + allResults[trialNum][n/2-1])/2;
   if(n%4 == 2){
     q1 = allResults[trialNum][n/4];
     q3 = allResults[trialNum][(3*n)/4];
   }else{
     q1 = (allResults[trialNum][n/4]+allResults[trialNum][n/4-1])/2;
     q3 = (allResults[trialNum][(3*n)/4]+allResults[trialNum][(3*n)/4-1])/2;
   }
  }else{
   median = allResults[trialNum][n/2]/1000;
   if(n%4 == 1){
     q1 = (allResults[trialNum][n/4]+allResults[trialNum][n/4-1])/2;
     q3 = (allResults[trialNum][(3*n)/4]+allResults[trialNum][(3*n)/4+1])/2;
   }else{
     q1 = allResults[trialNum][n/4];
     q3 = allResults[trialNum][(3*n)/4];
   }
  }
  
  cout << "Whisker plot (min-q1-median-q2-max): " << min << "--"<< q1 << "--" << median
   	<< "--" << q3 << "--" << max<< "(us)" << endl; 
  //cout << "Multiplied by num of batches: " << min*trialsToDo[trialNum]/1000 << "--"<<
 // 	 q1*trialsToDo[batchNum]/1000 << "--" << median*trialsToDo[trialNum]/1000 << "--"
 // 	 << q3*trialsToDo[batchNum]/1000 << "--" << max*trialsToDo[trialNum]/1000 <<
 // 	 "(ms)" << endl;
}

uint64_t sumAllTimes(int trialNum){
  int n = numUpdates*checks;
  uint64_t sum = 0;
  for(int i = 0; i< n; i++){
      sum += allResults[trialNum][i];
  }
  return sum;
}

double allResultsSTDEV(int trialNum){
  int n = numUpdates*checks;
  double sum = 0;
  for(int i = 0; i< n; i++){
      sum += allResults[trialNum][i];
  }

  double average = sum/n;
  double sumSqDiff = 0;
  
  for(int i = 0; i< n; i++){
      sumSqDiff += (allResults[trialNum][i] - average)*(allResults[trialNum][i] - average);
  }
  
  double ans = sqrt(sumSqDiff/(n - 1));
  
  return ans;
}

int getNonCheck(){
  int n;
  bool works = false;
  while(!works){
    n = rand();
    works = true;
    if(n > RAND_MAX - RAND_MAX%databaseSize) works = false;
    n = n%databaseSize;
    for(int i = 0; i< checks; i++){
      if(n == checkIdxs[i]) works = false;
    }
  }
  return n;
}

int getCheck(){
  return checkIdxs[rand()%checks];
}

int query_test_w_update(uint64_t num_items, uint64_t item_size, uint32_t degree,
               uint32_t lt, uint32_t dim) {
  uint64_t number_of_items = num_items;
  uint64_t size_per_item = item_size; // in bytes
  uint32_t N = degree;

  // Recommended values: (logt, d) = (12, 2) or (8, 1).
  uint32_t logt = lt;
  uint32_t d = dim;

  EncryptionParameters enc_params(scheme_type::bfv);
  PirParams pir_params;

  // Generates all parameters

  //cout << "Main: Generating SEAL parameters" << endl;
  gen_encryption_params(N, logt, enc_params);

  //cout << "Main: Verifying SEAL parameters" << endl;
  verify_encryption_params(enc_params);
  //cout << "Main: SEAL parameters are good" << endl;

  //cout << "Main: Generating PIR parameters" << endl;
  gen_pir_params(number_of_items, size_per_item, d, enc_params, pir_params);

  // gen_params(number_of_items, size_per_item, N, logt, d, enc_params,
  // pir_params);
  //print_pir_params(pir_params);

  //cout << "Main: Initializing the database (this may take some time) ..."
      // << endl;

  // Create test database
  auto db(make_unique<uint8_t[]>(number_of_items * size_per_item));
  //cout << "Database size: " << (db.get()).size();

  // Copy of the database. We use this at the end to make sure we retrieved
  // the correct element.
  auto db_copy(make_unique<uint8_t[]>(number_of_items * size_per_item));

  random_device rd;
  for (uint64_t i = 0; i < number_of_items; i++) {
    for (uint64_t j = 0; j < size_per_item; j++) {
      uint8_t val = rd() % 256;
      db.get()[(i * size_per_item) + j] = val;
      db_copy.get()[(i * size_per_item) + j] = val;
    }
  }

  // Initialize PIR Server
  //cout << "Main: Initializing server and client" << endl;
  PIRServer server(enc_params, pir_params);

  // Initialize PIR client....
  PIRClient client(enc_params, pir_params);
  GaloisKeys galois_keys = client.generate_galois_keys();

  // Set galois key for client with id 0
  //cout << "Main: Setting Galois keys...";
  server.set_galois_key(0, galois_keys);

  // Measure database setup
  auto time_pre_s = high_resolution_clock::now();
  server.set_database(move(db), number_of_items, size_per_item);
  server.preprocess_database();
  //cout << "Main: database pre processed " << endl;
  auto time_pre_e = high_resolution_clock::now();
  auto time_pre_us =
      duration_cast<microseconds>(time_pre_e - time_pre_s).count();
  
  auto total_time_query_us = 0;
  auto total_time_server_us = 0;
  auto total_time_decode_us = 0;
  //Create a query for every single check
  uint64_t indexes[checks];
  uint64_t offsets[checks];
  PirQuery queries[checks];
  
  for(int i = 0; i< checks; i++){
    indexes[i] = client.get_fv_index(checkIdxs[i]);
    offsets[i] = client.get_fv_offset(checkIdxs[i]);
    queries[i] = client.generate_query(indexes[i]);
    cout <<"Query size: " <<  queries[i].size() << endl;
  }
  
  //Old code for a single index:
  //uint64_t ele_index =
  //rd() % number_of_items; // element in DB at random position
  //uint64_t index = client.get_fv_index(ele_index);   // index of FV plaintext
  //uint64_t offset = client.get_fv_offset(ele_index); // offset in FV plaintext
  //cout << "Main: element index = " << ele_index << " from [0, "
  //  << number_of_items - 1 << "]" << endl;
  //cout << "Main: FV index = " << index << ", FV offset = " << offset << endl;
  
  // Measure query generation
  //auto time_query_s = high_resolution_clock::now();
  //PirQuery query = client.generate_query(index);
  //auto time_query_e = high_resolution_clock::now();
  //auto time_query_us =
  //duration_cast<microseconds>(time_query_e - time_query_s).count();
  //cout << "Main: query generated" << endl; 
  
  // To marshall query to send over the network, you can use
  // serialize/deserialize: std::string query_ser = serialize_query(query);
  // PirQuery query2 = deserialize_query(d, 1, query_ser, CIPHER_SIZE);
  
  //get and check initial replies
  PirReply replies[checks];
  vector<uint8_t> dec_replies[checks];
  
  for(int j = 0; j< checks; j++){
    replies[j] = server.generate_reply(queries[j],0);
    cout << "Reply size: " << replies[j].size() << endl;
    dec_replies[j] = client.decode_reply(replies[j], offsets[j]);
    cout << "Decoded Reply size: " << replies[j].size() << endl;
    assert(dec_replies[j].size() == size_per_item);
    
    bool failed = false;
    // Check that we retrieved the correct element
    for (uint32_t i = 0; i < size_per_item; i++) {
      if (dec_replies[j][i] != db_copy.get()[(checkIdxs[j] * size_per_item) + i]) {
        failed = true;
      }
    }
    if (failed) {
      cout << "Failed during inintial query at loop # " << j <<endl;
      return -1;
    }
  }
  
  for(int c = 0; c < numUpdates; c++){
    uint64_t changeIdx = 0;
    if(rand()%2 == 0){
      cout << "Changing check idx" << endl;
      changeIdx = getCheck();
    }else{
      cout << "Changing non-check" <<endl;
      changeIdx = getNonCheck();
    }
    cout << "Change idx: " << changeIdx <<  endl;
    
    auto newDbEntry(make_unique<uint8_t[]>(size_per_item));
    vector<uint8_t> text;
    for (uint64_t j = 0; j < size_per_item; j++) {
      uint8_t val = rd() % 256;
      newDbEntry.get()[j] = val;
      // changes database copy
      db_copy.get()[(changeIdx * size_per_item) + j] = val;
    }
    //need to change this to use text vector instead of a random number
    
    Plaintext pt = server.changed_db_at_idx(move(newDbEntry),changeIdx,size_per_item);
    for(int i = 0; i< checks; i++){
      Plaintext ptCopy = pt;
      server.update(queries[i], replies[i], changeIdx, ptCopy,0);
      dec_replies[i] = client.decode_reply(replies[i], offsets[i]);
    }
    
    server.simple_set(changeIdx/pir_params.elements_per_plaintext,pt);
    
    for(int j = 0; j< checks; j++){
      bool failed = false;
      
      if(dec_replies[j].size() != size_per_item){
        failed = true;
        cout << "Wrong size" << endl;
      }
      // Check that we retrieved the correct element
      for (uint32_t i = 0; i < size_per_item; i++) {
        if (dec_replies[j][i] != db_copy.get()[(checkIdxs[j] * size_per_item) + i]) {
          failed = true;
        }
      }
      
      if(failed){
        cout << "Failed on the " << c<< "th update, changed the " << changeIdx<< "th index, and failed on the " << j <<"th check." << endl;
        return -1;
      }
    }
  }
 
  
  // Measure query processing (including expansion)
  //auto time_server_s = high_resolution_clock::now();
  //PirReply reply = server.generate_reply(query, 0);
  //auto time_server_e = high_resolution_clock::now();
  //auto time_server_us =bool correct_elem(vector<uint8_t> dec_reply, int index)
  //duration_cast<microseconds>(time_server_e - time_server_s).count();
  
  // Measure response extraction
  //auto time_decode_s = chrono::high_resolution_clock::now();
  //vector<uint8_t> elems = client.decode_reply(reply, offset);
  //auto time_decode_e = chrono::high_resolution_clock::now();
  //auto time_decode_us =
  //duration_cast<microseconds>(time_decode_e - time_decode_s).count();
  
  //assert(elems.size() == size_per_item);

  //bool failed = false;
  // Check that we retrieved the correct element
  //for (uint32_t i = 0; i < size_per_item; i++) {
  //  if (elems[i] != db_copy.get()[(ele_index * size_per_item) + i]) {
      //cout << "Main: elems " << (int)elems[i] << ", db "
      //     << (int)db_copy.get()[(ele_index * size_per_item) + i] << endl;
      //cout << "Main: PIR result wrong at " << i << endl;
  //    failed = true;
  //  }
  //}
  //if (failed) {
  //  return -1;
 // }
  
  //--------------------------Updates and re-checking:---------------------------

  // add times to total time
  //total_time_query_us += time_query_us;
  //total_time_server_us += time_server_us;
  //total_time_decode_us += time_decode_us;

  // Output results
  //cout << "Main: PIR result correct!" << endl;
  //cout << "Main: PIRServer pre-processing time: " << time_pre_us / 1000 << " ms"
  //     << endl;
  //cout << "Main: PIRClient query generation time: " << total_time_query_us/
  //	trialsPerDatabase / 1000 << " ms" << endl;
  //cout << "Main: PIRServer reply generation time: " << total_time_server_us/
  //	trialsPerDatabase / 1000 << " ms" << endl;
  //cout << "Main: PIRClient answer decode time: " << total_time_decode_us/
  //	trialsPerDatabase / 1000<< " ms" << endl;
  //cout << "Main: Reply num ciphertexts: " << reply.size() << endl;
  return 0;
}
