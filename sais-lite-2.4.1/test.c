/*
 * test.c for sais-lite
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sais.h"


int
main(int argc, const char *argv[]) {

  char text[] = "banana";
  int text_len;
  int *suffix_array;
  int *inverse_sa;

  text_len = strlen(text);

  suffix_array = malloc(text_len * sizeof(int));
  inverse_sa = malloc(text_len * sizeof(int));


  if (sais((unsigned char *) text, suffix_array, text_len) != 0) {
    printf("sais failed\n");
    exit(EXIT_FAILURE);
  }

  printf("text:\n\t %s\n\n", text);
  
  printf("SA indices: {\n\t");
  for (int i = 0; i < text_len; i++) {
    printf("%d, ", suffix_array[i]);
  }
  printf("\n}\n");

  for (int i=0; i < text_len; i++) {
    inverse_sa[suffix_array[i]] = i;
  }

  printf("ISA indices: {\n\t");
  for (int i = 0; i < text_len; i++) {
    printf("%d, ", inverse_sa[i]);
  }
  printf("\n}\n");



  return 0;
}
