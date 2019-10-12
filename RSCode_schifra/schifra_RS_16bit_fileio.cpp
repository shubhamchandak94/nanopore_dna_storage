/*
(**************************************************************************)
(*                                                                        *)
(*                                Schifra                                 *)
(*                Reed-Solomon Error Correcting Code Library              *)
(*                                                                        *)
(* Release Version 0.0.1                                                  *)
(* http://www.schifra.com                                                 *)
(* Copyright (c) 2000-2019 Arash Partow, All Rights Reserved.             *)
(*                                                                        *)
(* The Schifra Reed-Solomon error correcting code library and all its     *)
(* components are supplied under the terms of the General Schifra License *)
(* agreement. The contents of the Schifra Reed-Solomon error correcting   *)
(* code library and all its components may not be copied or disclosed     *)
(* except in accordance with the terms of that agreement.                 *)
(*                                                                        *)
(* URL: http://www.schifra.com/license.html                               *)
(*                                                                        *)
(**************************************************************************)
*/


/*
   Description: This example will demonstrate how to instantiate a 16-bit per
                symbol Reed-Solomon encoder and decoder, add the full amount
                of possible errors, correct the errors, and output the various
                pieces of relevant information. Furthermore this example will
                demonstrate the use of the Reed-Solomon codec without the use
                of LUTs for the finite field computations.
*/


#include <cstddef>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#define NO_GFLUT
#include "schifra_galois_field.hpp"
#undef NO_GFLUT
#include "schifra_galois_field_polynomial.hpp"
#include "schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra_reed_solomon_encoder.hpp"
#include "schifra_reed_solomon_decoder.hpp"
#include "schifra_reed_solomon_block.hpp"
#include "schifra_error_processes.hpp"

#include "schifra_reed_solomon_bitio.hpp"

#include "RS_paramaters_from_python.hpp"



int main(int argc, char* argv[])
{
   /* if argv[1] = 1, we encode input_data. if argv[1] = 0, we decode input_data  */


   /* Instantiate Finite Field and Generator Polynomials */
   const schifra::galois::field field(field_descriptor,
                                      schifra::galois::primitive_polynomial_size14,
                                      schifra::galois::primitive_polynomial14);

   schifra::galois::field_polynomial generator_polynomial(field);

   if (
        !schifra::make_sequential_root_generator_polynomial(field,
                                                            generator_polynomial_index,
                                                            generator_polynomial_root_count,
                                                            generator_polynomial)
      )
   {
      std::cout << "Error - Failed to create sequential root generator!" << std::endl;
      return 1;
   }

	
   /* if argv[1] = true, we encode input_data. if argv[1] = false, we decode input_data  */
   bool enc;
   std::istringstream(argv[1]) >> enc;
   
   // taking in file names
   const std::string input_file_name     = argv[2];
   const std::string output_file_name = argv[3];
   
   std::cout << input_file_name << "\n" << output_file_name << "\n";
   
      
   // instnatiate a schifra block (which is what we do encoding and decoding on)
   schifra::reed_solomon::block<code_length,fec_length> block;   
   
   // reading input file via binary io and storing as a vec
   std::ifstream myinputfile (input_file_name, std::ios::binary);

   std::vector<uint16_t> messagevec;
		
   uint16_t tempmsg;
   myinputfile.read((char*)&tempmsg, sizeof(uint16_t));
   while(!myinputfile.eof()) 
   {
	  messagevec.push_back(tempmsg);
	  myinputfile.read((char*)&tempmsg, sizeof(uint16_t));
   }
   myinputfile.close();


   // converting vector to "block" data type which is what schifra uses
   schifra::reed_solomon::copy(&messagevec[0], messagevec.size(), block);
      
      
   int output_length; 
   
   if(enc) /* this means we want to encode */
   {
	   /* Instantiate Encoder */
	  typedef schifra::reed_solomon::encoder<code_length,fec_length,data_length> encoder_t;
	  const encoder_t encoder(field, generator_polynomial);

	  output_length = code_length;

	   /* Transform block into Reed-Solomon encoded codeword */
	   if (!encoder.encode(block))
	   {
		  std::cout << "Error - Critical encoding failure! "
	     			<< "Msg: " << block.error_as_string()  << std::endl;
	      return 1;
	   }

	  

	   
	   std::cout << "Encoder Parameters [" << encoder_t::trait::code_length << ","
	    								   << encoder_t::trait::data_length << ","
		    							   << encoder_t::trait::fec_length  << "]" << std::endl;
	
   }



   
   if(!enc) /* this means we want to decode */
   { 
	  /* Instantiate Decoder */
	  typedef schifra::reed_solomon::decoder<code_length,fec_length,data_length> decoder_t;
	  const decoder_t decoder(field, generator_polynomial_index);

	  output_length = data_length;

	  // deciding whether or not we have erasures
	  /* if argv[4] = true, we have erasures. if argv[4] = false, we do not have erasures  */
      bool erasure_true;
      std::istringstream(argv[4]) >> erasure_true;

      // taking in erasure_lication_list file name   
      const std::string erasure_file_name  = argv[5];   

	  // reading in erasure list
	  std::vector<size_t> erasure_location_list;
	
	  // if we do have erasures do the following
      if(erasure_true)
      {
	     std::ifstream erasureinputfile(erasure_file_name,std::ios::binary);

	     uint16_t temp_erasure_loc;
	     erasureinputfile.read((char*)&temp_erasure_loc, sizeof(uint16_t));
	     while(!erasureinputfile.eof()) 
	     {
		    erasure_location_list.push_back((size_t)temp_erasure_loc);
		    erasureinputfile.read((char*)&temp_erasure_loc, sizeof(uint16_t));
	     }
	     erasureinputfile.close();
     

      if (!decoder.decode(block,erasure_location_list))
      {
          std::cout << "Error - Critical decoding failure!"
                    << "Msg: " << block.error_as_string()  << std::endl;
          return 1;
      }
	
	}


	// if no erasures, do the following
	if(!erasure_true)
	{
     if (!decoder.decode(block))
     {
      std::cout << "Error - Critical decoding failure!"
                << "Msg: " << block.error_as_string()  << std::endl;
      return 1;
     }
    }
	  
	  
	std::cout << "Decoder Parameters [" << decoder_t::trait::code_length << ","
	     								<< decoder_t::trait::data_length << ","
										<< decoder_t::trait::fec_length  << "]" << std::endl;
      
   }

	 // writing output to a file via binary io
	  std::ofstream myoutputfile(output_file_name,std::ios::binary);
	  uint16_t  tmpthing;
  
	  for ( int i = 0; i < output_length; i++) 
	  {
		tmpthing = block.data[i]&0xFFFF;
		myoutputfile.write((char*)&tmpthing,sizeof(uint16_t));
	  }


	  myoutputfile.close();


   return 0;
}




