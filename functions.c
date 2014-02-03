#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "functions.h"
char * read_is_element(char * seqname)
// Ielasa IS6110.txt fasta formātā un atgriež dns sekvenci
{ 
  char isname[1000]="";
  strcat(isname, seqname);
  static char * sequence;
  FILE *filepointer; 
  //static char sequence, *temp;
  static char *temp;
  sequence=malloc(sizeof(char));
  if(sequence==NULL){
    printf("ERROR allocating memory for sequence\n");
  }
  int is_sequence=0;
  int character;
  int nnuc=0;
  filepointer=fopen(isname,"read"); //"/home/ivars/Learn/C/find_IS/IS6110.txt"
  while ((character=fgetc(filepointer))!=EOF){
    if (is_sequence==1 && character!='\n' && character!=EOF){
      sequence[nnuc]=character;
      temp=realloc(sequence, (nnuc+2)*sizeof(char));
      if (temp !=NULL){
        sequence=temp;
      } else {
        free(sequence);
        printf("Error allocating memory!\n");
        //return 1;
      }
      ++nnuc;
    }
    if (character=='>'){
      is_sequence=0;
    }
    if (character=='\n' && is_sequence==0){
      is_sequence=1;
    }
  }
  sequence[nnuc]='\0';
  fclose(filepointer); 
  return sequence;
}

void get_fragment(int start, int end,char * fragment, char * sequence)  //aizvietots ar library funkciju strncpy
// funkcija kas atgriež fragment no * sequence pointera, ja doti arraya (pointer) indexi
{  
  //static char * fragment;
  int i=0;
  //fragment=malloc((end-start+1)*sizeof(char));
  if (fragment==NULL){
    printf("Error allocating memory at get_fragment\n");
  }
  //printf("%p\n",fragment);
  memset(fragment,'\0',end-start+1);
  for(start;start<=end;++start){
    fragment[i]=sequence[start];
    ++i;
  }
  fragment[i]='\0';
  
 // return fragment;
}

void create_revcom(char * revcom,char * sequence, int last_index)
// izveido sekvencei reverse complement
{
  int i;
  int j=0;
//  static char * revcom;
//  revcom = malloc((last_index+2)*sizeof(char));
  for (i=last_index;i>=0;--i){
    if (sequence[i]=='A'){
      revcom[j]='T';
    } else if  (sequence[i]=='T'){
      revcom[j]='A';
    } else if (sequence[i]=='G'){
      revcom[j]='C';
    } else if (sequence[i]=='C'){
      revcom[j]='G';
    } else if (sequence[i]=='N'){
      revcom[j]='N';
    } else {
      printf("%c nav standarta nukleotīds!\n",revcom[j]);
    }
    ++j;
  }
  revcom[j]='\0';
  //return revcom;
}
// reverse the given null-terminated string in place
void inplace_reverse(char * str)
{
  if (str)
  {
    char * end = str + strlen(str) - 1;

    // swap the values in the two given variables
    // XXX: fails when a and b refer to same memory location
#   define XOR_SWAP(a,b) do\
    {\
      a ^= b;\
      b ^= a;\
      a ^= b;\
    } while (0)

    // walk inwards from both ends of the string, 
    // swapping until we get to the middle
    while (str < end)
    {
      XOR_SWAP(*str, *end);
      str++;
      end--;
    }
#   undef XOR_SWAP
  }
}
/*
char * extract_frag(char *sequence, int is_start, int coord){
//  char * sequence =malloc(1+strlen(seq));
//  if (sequence){
 //   sequence=strdup ( seq);
//  } else {
//    printf ("malloc failure!\n");
//  }
  int i;
  int j=0;
  int length=0;
  static char * frag, * temp;
  frag=malloc(2*sizeof(char));
  while (sequence[length]!='\0'){
    ++length;
  }
  if (is_start==1){
    for(i=coord+1;i<=length;++i){
      frag[j]=sequence[i];
      //allocate memory to frag
      temp=realloc(frag, (j+2)*sizeof(char));
      if (temp!=NULL){
        frag=temp;
      } else {
        free(frag);
        printf("Error allocating memory for input string frag!\n");
      }
      // done
    ++j;
    }
    frag[j]='\0';
    return frag;
  
  } else if (is_start==0){
    for(i=0;i<coord;++i){
      frag[i]=sequence[i];
          //allocate memory to frag
      temp=realloc(frag, (j+2)*sizeof(char));
      if (temp!=NULL){
        frag=temp;
      } else {
        free(frag);
        printf("Error allocating memory for input string frag!\n");
      }
    }
    // done
    frag[i]='\0';
    return frag;
    
  }
}*/
