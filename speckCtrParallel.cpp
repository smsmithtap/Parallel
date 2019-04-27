/*
 * Copyright (c) 2016-2017 Naruto TAKAHASHI <tnaruto@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <random>
#include "speck.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

// https://eprint.iacr.org/2013/404.pdf
//
// Speck128/256
//  Key:        1f1e1d1c1b1a1918 1716151413121110 0f0e0d0c0b0a0908 0706050403020100
//  Plaintext:  65736f6874206e49 202e72656e6f6f70
//  Ciphertext: 4109010405c0f53e 4eeeb48d9c188f43
static const uint64_t s_plain_text[2] = {0x202e72656e6f6f70, 0x65736f6874206e49};
static const uint64_t s_cipher_text[2] = {0x4eeeb48d9c188f43, 0x4109010405c0f53e};

static const uint8_t s_key_stream[32] = {
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
};

//Set to "The test string/".  plain_text_stream Will be concatenated together Block_num times each loop in the non random stream test.  So twice on the second, 
//3 times on the third, etc.
static const uint8_t s_plain_text_stream[16] = {
	0x54, 0x68, 0x65, 0x20, 0x74, 0x65, 0x73, 0x74, 0x20, 0x73, 0x74, 0x72, 0x69, 0x6e, 0x67, 0x2f
};
static const uint8_t s_cipher_text_stream[16] = {
    0x43, 0x8f, 0x18, 0x9c, 0x8d, 0xb4, 0xee, 0x4e, 0x3e, 0xf5, 0xc0, 0x05, 0x04, 0x01, 0x09, 0x41,
};

//TEST
static const uint8_t s_nonce[16] = {
    0x45, 0x8f, 0x18, 0x9c, 0x8d, 0xb4, 0xee, 0x4e, 0x3e, 0xf5, 0xc0, 0x05, 0x04, 0x01, 0x09, 0x41,
};

int thread_count;

#define BLOCK_SIZE 16

FILE* OUTPUTFILE;

void generate_random_array(uint8_t *iv, size_t iv_len) {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis;

    for(int i=0; i<iv_len; i++) {
        iv[i] = static_cast<uint8_t>(dis(gen));
    }
}

void show_array(const char *explain, const uint8_t *array, size_t len) {
    printf("%20s ", explain);
    for(int i=len-1; i >= 0; i--) {
        printf("%02x ", array[i]);
    }
    printf("\n");
}

//Output ASCII Hex values to an output file, to show the plain text and decrypted stream are the same.
void output_array_to_file(const char *explain, const uint8_t *array, size_t len, FILE* output)
{
	fprintf(output, "%20s ", explain);
    for(int i=len-1; i >= 0; i--) {
        fprintf(output, "%02x ", array[i]);
    }
    fprintf(output, "\n");
}

int encrypt_decrypt_stream_test(int block_num) {
    int r = 0;
    speck_ctx_t *ctx = NULL;
    uint8_t *plain_text_stream = NULL;
    uint8_t *crypted_text_stream = NULL;
    uint8_t *decrypted_text_stream = NULL;
    uint8_t *iv_text_stream = NULL;
    uint8_t *origin_iv_text_stream = NULL;
	
	//Calculate chunk size and remainder for splitting plain text blocks between threads.
	int blockChunk = (block_num/thread_count);
	int blockChunkRemainder = (block_num%thread_count);

    plain_text_stream = (uint8_t*)malloc(BLOCK_SIZE * block_num);
    if (!plain_text_stream) {
        r = 1;
        goto finish;
    }
    crypted_text_stream = (uint8_t*)malloc(BLOCK_SIZE * block_num);
    if (!crypted_text_stream) {
        r = 1;
        goto finish;
    }
    decrypted_text_stream = (uint8_t*)malloc(BLOCK_SIZE * block_num);
    if (!decrypted_text_stream) {
        r = 1;
        goto finish;
    }
    iv_text_stream = (uint8_t*)malloc(BLOCK_SIZE);
    if (!iv_text_stream) {
        r = 1;
        goto finish;
    }
    origin_iv_text_stream = (uint8_t*)malloc(BLOCK_SIZE);
    if (!origin_iv_text_stream) {
        r = 1;
        goto finish;
    }
	
	//Set up nonce (random number that is not reused) in origin_iv_text_stream to serve as additional key for counter mode block cipher.
	//generate_random_array(origin_iv_text_stream, BLOCK_SIZE);
	
	//TEST - Keep nonce the same for testing.
	memcpy(origin_iv_text_stream, s_nonce, BLOCK_SIZE);

    for (int i = 0; i < block_num; i++) {
        memcpy(plain_text_stream + (i * BLOCK_SIZE), s_plain_text_stream, sizeof(s_plain_text_stream));
    }
	
	//show_array("plain", plain_text_stream, block_num * BLOCK_SIZE);
	
	fprintf(OUTPUTFILE, "Plain Text Before Encryption:  %s\n", plain_text_stream);
	//Output ASCII Hex values to an output file, to show the plain text and decrypted stream are the same.
	output_array_to_file("Plain Text Hex:", plain_text_stream, block_num*BLOCK_SIZE, OUTPUTFILE);
	printf("Size of plain_text_stream:  %d.\n", block_num*BLOCK_SIZE);
	
	printf("Block num:  %d.  Size of s_plain_text_stream:  %d.\n", block_num, sizeof(s_plain_text_stream));
	printf("Size of plain_text_stream:  %d.  Size of uint8_t:  %d.  Length of plain_text_stream:  %d.\n", sizeof(plain_text_stream), sizeof(uint8_t), sizeof(plain_text_stream)/sizeof(uint8_t));
	
	fflush(stdout);
	
	//Moved speck_init above parallel section so each process can use the same ctx value.
	ctx = speck_init(SPECK_ENCRYPT_TYPE_128_256, s_key_stream, sizeof(s_key_stream));
	
	if (!ctx) {
        r = 1;
        goto finish;
    }

	//Begin parallel section for encryption.
	//First determine number of blocks in the plain_text_stream, then split the blocks out to each process.
	//16 Blocks = 64 bits with the uint8_t data type.
	//Shared memory with OpenMP, so no need to use separate variables.
	#pragma omp parallel num_threads(thread_count)
	{
		int cur_thread,tot_thread_count;
		//Retrieve the current thread ID.
		cur_thread = omp_get_thread_num();
		//Retrieve the total number of threads.
		tot_thread_count = omp_get_num_threads();
		
		int first = (blockChunk*cur_thread)*BLOCK_SIZE;
		int curBlockChunkSize = blockChunk;
		
		if(cur_thread == (tot_thread_count-1))
		{
			curBlockChunkSize += blockChunkRemainder;
		}
		
		//if(first > 0)
		//{
		//	first = first-1;
		//}
		
		//May need to move this out?  Or may need to keep it for each process.
		//ctx = speck_init(SPECK_ENCRYPT_TYPE_128_256, s_key_stream, sizeof(s_key_stream));
		//if (!ctx) {
		//    r = 1;
		//    goto finish;
		//}
		
		memcpy(iv_text_stream, origin_iv_text_stream, BLOCK_SIZE);
		
		fprintf(OUTPUTFILE, "Nonce:  %s\n", iv_text_stream);
		//Serial version of the encryption call.
		//r = speck_ctr_encrypt(ctx, plain_text_stream, crypted_text_stream, BLOCK_SIZE * block_num, iv_text_stream, BLOCK_SIZE);
		//if(cur_thread == 1)
		//{
			printf("Calling encryption function - thread num:  %d first index in text stream:  %d.  Length of plain text to encrypt: %d\n", cur_thread, first, int(BLOCK_SIZE * curBlockChunkSize));
			fflush(stdout);
		
			printf("Thread %d calling encryption with Plain_text_stream address:  %p. Crypted_text_stream address:  %p.\n", cur_thread, (plain_text_stream+(first)), (crypted_text_stream+(first)));
			r = speck_ctr_encrypt(ctx, (plain_text_stream+(first)), (crypted_text_stream+(first)), BLOCK_SIZE * curBlockChunkSize, iv_text_stream, BLOCK_SIZE, cur_thread);
		//}
		
		printf("Finished with encryption function - thread num:  %d first index in text stream:  %d.\n", cur_thread, first);
		fflush(stdout);
	}
	//End parallel encryption section.
	
    if (r < 0) {
        r = 1;
        goto finish;
    }
	
	fprintf(OUTPUTFILE, "Nonce after parallel:  %s\nEncrypted text:  %s\n", origin_iv_text_stream, crypted_text_stream);
	//output_array_to_file("Nonce Hex Values after parallel:", iv_text_stream, BLOCK_SIZE, OUTPUTFILE);
	
	fflush(OUTPUTFILE);
	
	output_array_to_file("Encrypted Text Hex Values:", crypted_text_stream, block_num*BLOCK_SIZE, OUTPUTFILE);
	fflush(OUTPUTFILE);
	
	//Begin parallel section for decryption.
	//First determine number of blocks in the plain_text_stream, then split the blocks out to each process.
	//Shared memory with OpenMP, so no need to use separate variables.
	/*
	#pragma omp parallel num_threads(thread_count)
	{
		int cur_thread,tot_thread_count;
		//Retrieve the current thread ID.
		cur_thread = omp_get_thread_num();
		//Retrieve the total number of threads.
		tot_thread_count = omp_get_num_threads();
		
		int first = (blockChunk*cur_thread)*BLOCK_SIZE;
		int curBlockChunkSize = blockChunk;
		
		if(cur_thread == (tot_thread_count-1))
		{
			curBlockChunkSize += blockChunkRemainder;
		}
		
		//if(first > 0)
		//{
		//	first = first-1;
		//}
		
		//May need to move this out?  Or may need to keep it for each process.
		//ctx = speck_init(SPECK_ENCRYPT_TYPE_128_256, s_key_stream, sizeof(s_key_stream));
		//if (!ctx) {
		//    r = 1;
		//    goto finish;
		//}
		
		memcpy(iv_text_stream, origin_iv_text_stream, BLOCK_SIZE);
		
		fprintf(OUTPUTFILE, "Nonce:  %s\n", iv_text_stream);
		//Serial version of the encryption call.
		//r = speck_ctr_encrypt(ctx, plain_text_stream, crypted_text_stream, BLOCK_SIZE * block_num, iv_text_stream, BLOCK_SIZE);
		printf("Calling decryption function - thread num:  %d first index in text stream:  %d.\n", cur_thread, first);
		fflush(stdout);
		
		r = speck_ctr_encrypt(ctx, crypted_text_stream+(first*BLOCK_SIZE), decrypted_text_stream+(first*BLOCK_SIZE), BLOCK_SIZE * curBlockChunkSize, iv_text_stream, BLOCK_SIZE);
		
		printf("Finished with decryption function - thread num:  %d first index in text stream:  %d.\n", cur_thread, first);
		fflush(stdout);
	}
	*/
	
    memcpy(iv_text_stream, origin_iv_text_stream, BLOCK_SIZE);
    r = speck_ctr_decrypt(ctx, crypted_text_stream, decrypted_text_stream, BLOCK_SIZE * block_num, iv_text_stream, BLOCK_SIZE, 0);
	
	//End parallel decryption section.
	
    if (r < 0) {
        r = 1;
        goto finish;
    }
	
	fprintf(OUTPUTFILE, "Nonce after decryption:  %s\nEncrypted Text After Decryption:  %s\n", origin_iv_text_stream, crypted_text_stream);
	
	output_array_to_file("Encrypted Text Hex Values After Decryption:", crypted_text_stream, block_num*BLOCK_SIZE, OUTPUTFILE);
	fflush(OUTPUTFILE);
	
	//Output ASCII Hex values to an output file, to show the plain text and decrypted stream are the same.
	output_array_to_file("Plain Text Hex:", plain_text_stream, block_num*BLOCK_SIZE, OUTPUTFILE);
	
	output_array_to_file("Decrypted text Hex:", decrypted_text_stream, block_num*BLOCK_SIZE, OUTPUTFILE);
	
	fprintf(OUTPUTFILE, "Plain Text Before Encryption:  %s\n", plain_text_stream);
	fprintf(OUTPUTFILE, "Decrypted text:  %s\n\n", decrypted_text_stream);
	
	fflush(OUTPUTFILE);
	
	//Testing to ensure the plain_text_stream and decrypted_text_stream are the same.  This shows that the decryption was completed successfully.
    for (int i = 0; i < BLOCK_SIZE * block_num; i++) {
        if (plain_text_stream[i] != decrypted_text_stream[i]) {
            printf("block_num:%d idx:%d  0x%02x != 0x%02x\n", block_num, i, plain_text_stream[i], decrypted_text_stream[i]);
            show_array("iv", origin_iv_text_stream, BLOCK_SIZE);
            show_array("plain", plain_text_stream, block_num * BLOCK_SIZE);
            show_array("decrypted", decrypted_text_stream, block_num * BLOCK_SIZE);
            show_array("counted iv", iv_text_stream, BLOCK_SIZE);
            printf("\n");

            r = 1;
            goto finish;
        }
    }

    finish:
    free(plain_text_stream);
    free(crypted_text_stream);
    free(decrypted_text_stream);
    free(iv_text_stream);
    free(origin_iv_text_stream);

    speck_finish(&ctx);
    return r;
}

int encrypt_decrypt_random_stream_test(int block_num) {
    int r = 0;
    speck_ctx_t *ctx = NULL;
    uint8_t *key_text_stream = NULL;
    int key_text_length = sizeof(s_key_stream);
    uint8_t *plain_text_stream = NULL;
    uint8_t *crypted_text_stream = NULL;
    uint8_t *decrypted_text_stream = NULL;
    uint8_t *iv_text_stream = NULL;
    uint8_t *origin_iv_text_stream = NULL;

    key_text_stream = (uint8_t*)malloc(key_text_length);
    if(!key_text_stream) {
        r = 1;
        goto finish;
    }
    plain_text_stream = (uint8_t*)malloc(block_num);
    if (!plain_text_stream) {
        r = 1;
        goto finish;
    }
    crypted_text_stream = (uint8_t*)malloc(block_num);
    if (!crypted_text_stream) {
        r = 1;
        goto finish;
    }
    decrypted_text_stream = (uint8_t*)malloc(block_num);
    if (!decrypted_text_stream) {
        r = 1;
        goto finish;
    }
    iv_text_stream = (uint8_t*)malloc(BLOCK_SIZE);
    if (!iv_text_stream) {
        r = 1;
        goto finish;
    }
    origin_iv_text_stream = (uint8_t*)malloc(BLOCK_SIZE);
    if (!origin_iv_text_stream) {
        r = 1;
        goto finish;
    }

    generate_random_array(key_text_stream, key_text_length);
    generate_random_array(plain_text_stream, block_num);
    generate_random_array(origin_iv_text_stream, BLOCK_SIZE);

    ctx = speck_init(SPECK_ENCRYPT_TYPE_128_256, key_text_stream, key_text_length);
    if (!ctx) {
        r = 1;
        goto finish;
    }
    memcpy(iv_text_stream, origin_iv_text_stream, BLOCK_SIZE);
    r = speck_ctr_encrypt(ctx, plain_text_stream, crypted_text_stream, block_num, iv_text_stream, BLOCK_SIZE, 0);
    if (r < 0) {
        r = 1;
        goto finish;
    }
	
	//output_array_to_file("Plain Text - Hex ", plain_text_stream, block_num*BLOCK_SIZE, OUTPUTFILE);
	
	fprintf(OUTPUTFILE, "Decrypted text:  %s\n", crypted_text_stream);
	
	fflush(OUTPUTFILE);
	
    memcpy(iv_text_stream, origin_iv_text_stream, BLOCK_SIZE);
    r = speck_ctr_decrypt(ctx, crypted_text_stream, decrypted_text_stream, block_num, iv_text_stream, BLOCK_SIZE, 0);
    if (r < 0) {
        r = 1;
        goto finish;
    }
	
	//Output text form of plain text and decrypted text to show they are the same.
	fprintf(OUTPUTFILE, "Plain Text Before Encryption:  %s\n", plain_text_stream);
	fprintf(OUTPUTFILE, "Decrypted text:  %s\n\n", decrypted_text_stream);
	
	//output_array_to_file("Decrypted text - Hex ", decrypted_text_stream, block_num*BLOCK_SIZE, OUTPUTFILE);
	
	fflush(OUTPUTFILE);
	
    for (int i = 0; i < block_num; i++) {
        if (plain_text_stream[i] != decrypted_text_stream[i]) {
            printf("block_num:%d idx:%d  0x%02x != 0x%02x\n", block_num, i, plain_text_stream[i], decrypted_text_stream[i]);
            show_array("iv", origin_iv_text_stream, BLOCK_SIZE);
            show_array("plain", plain_text_stream, block_num);
            show_array("decrypted", decrypted_text_stream, block_num);
            show_array("counted iv", iv_text_stream, BLOCK_SIZE);
            printf("\n");

            r = 1;
            goto finish;
        }
    }

    finish:
    free(key_text_stream);
    free(plain_text_stream);
    free(crypted_text_stream);
    free(decrypted_text_stream);
    free(iv_text_stream);
    free(origin_iv_text_stream);

    speck_finish(&ctx);
    return r;
}

int main(int args, char * argv[]) {
	
	double runTimeStart = omp_get_wtime();
	
	if(args < 2)
	{
			fprintf(stderr, "Usage: ./program <NumberOfThreads>\n");
			return -1;
	}
	
	thread_count = atoi(argv[1]);
	
	//int blockNum = (4*1024)/BLOCK_SIZE;
	
	int blockNum = 128/BLOCK_SIZE;
	
	OUTPUTFILE = fopen("output.txt", "w");
    printf("test encrypt_decrypt_stream_test\n");
    //for (int i = 1024; i <(4*1024)/BLOCK_SIZE; i++) {
	//	fprintf(OUTPUTFILE, "Loop %d\n", i);
        int r = encrypt_decrypt_stream_test(blockNum);
        if(r != 0) {
            return r;
        }
    //}
    printf("success encrypt_decrypt_stream_test\n");

    printf("test encrypt_decrypt_random_stream_test\n");
    //for (int i = 0; i <(4*1024); i++) {
	//	fprintf(OUTPUTFILE, "Loop %d\n", i);
        r = encrypt_decrypt_random_stream_test(blockNum);
        if(r != 0) {
            return r;
        }
    //}
    printf("success encrypt_decrypt_random_stream_test\n");

	fclose(OUTPUTFILE);
	
	printf("Run Time Start:  %lf\n", runTimeStart);

	double runTimeEnd = omp_get_wtime();

	printf("Run Time End:  %lf\n", runTimeEnd);

	double totalRunTime = runTimeEnd-runTimeStart;

	printf("Total Run Time:  %lf\n", totalRunTime);
	fflush(stdout);
	
    return 0;
}