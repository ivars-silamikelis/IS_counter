#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define BUF_SIZE 1024
char * read_is_element(char * seqname);
void get_fragment(int start, int end,char * fragment, char * sequence);
void create_revcom(char * revcom,char * sequence, int last_index);
void inplace_reverse(char * str);
int main(int argc, char *argv[])
{ 
  if (argc==1){
    printf("Please provide output file prefix and insertion sequence filename\n\n");
    printf("-----------EXAMPLE---------------\n\n");
    printf("./IS_find sample_name IS6110.txt\n");
    return 0;
  }
  char start_filename[200]="";
  char end_filename[200]="";
  char * is_name=strdup(argv[2]);
  //printf("%s\n", argv[1]);
  strcat(start_filename, argv[1]);
  //printf("%s\n", argv[1]);
  strcat(end_filename, argv[1]);
  //printf("%s\n", argv[1]);
  strcat(start_filename, ".starts.fa");
  strcat(end_filename, ".ends.fa");
  FILE * fs=fopen(start_filename,"w");
  FILE * fe=fopen(end_filename,"w");
  int i,j;
  int m=1;
  int n=1;
  int progress=0;
  int small_progress;
  char * sequence;
  char * start_fragment=malloc(sizeof(char));
  char * end_fragment=malloc(sizeof(char));
  char * start_fragment_revcom=malloc(sizeof(char));
  char * end_fragment_revcom=malloc(sizeof(char));
  char * read_sequence=malloc((BUF_SIZE+1)*sizeof(char));
  char * read_start=malloc(sizeof(char));
  char * read_end=malloc(sizeof(char));
  char * read_start_good=malloc(sizeof(char));
  char * read_end_good=malloc(sizeof(char));
  char * read_start_good_revcom=malloc(sizeof(char));
  char * read_end_good_revcom=malloc(sizeof(char));
  char * rev_rev_start=malloc(sizeof(char));
  char * rev_rev_end=malloc(sizeof(char));
  int fasta_last_index;
  int read_last_nuc_index; 


  //Atmiņā tiek ielasīta IS6110 sekvence
  sequence=read_is_element(is_name);
  fasta_last_index=strlen(sequence)-1;
  //Atmiņā tiek ielasītas sekvences no bam vai fastq faila
  while(fgets (read_sequence, BUF_SIZE, stdin)){
    //Tiek parādīts cik daudz sekvenču no bam vai fastq faila ir apstrādātas
    ++small_progress;
    if (small_progress==100000){
      progress+=small_progress;
      printf("%d reads processed\n", progress);
      small_progress=0;
    }
    read_last_nuc_index=strlen(read_sequence)-2;
    j=149;
    if (j>read_last_nuc_index){
      j=read_last_nuc_index-5; //Minimālais flankējošās sekvences garums 5
    }
    if (j<=0){
      continue;
    }
    for(i=j;i>=15;--i){
      // Iegūst fragmentus no IS6110 sekvences
      // Mazāk piešķirtās atmiņas izraisa segfault
      //Jāveic check for errors
      start_fragment=realloc(start_fragment,(i+10)*sizeof(char));
      end_fragment=realloc(end_fragment,(i+10)*sizeof(char));
      start_fragment_revcom=realloc(start_fragment_revcom,(i+10)*sizeof(char));
      end_fragment_revcom=realloc(end_fragment_revcom,(i+10)*sizeof(char));
      read_start=realloc(read_start,(i+10)*sizeof(char));
      read_end=realloc(read_end,(i+10)*sizeof(char));

      read_start_good=realloc(read_start_good,(read_last_nuc_index+10)*sizeof(char));
      read_end_good=realloc(read_end_good,(read_last_nuc_index+10)*sizeof(char));
      read_start_good_revcom=realloc(read_start_good_revcom,(read_last_nuc_index+10)*sizeof(char));
      read_end_good_revcom=realloc(read_end_good_revcom,(read_last_nuc_index+10)*sizeof(char));
      rev_rev_start=realloc(rev_rev_start,(read_last_nuc_index+10)*sizeof(char));
      rev_rev_end=realloc(rev_rev_end,(read_last_nuc_index+10)*sizeof(char));

    // Revcom jāizveido 
      

      strncpy(start_fragment, sequence,i+1);
      start_fragment[i+1]='\0';
      strncpy(end_fragment, sequence+fasta_last_index-i,i+1);
      end_fragment[i+1]='\0';
 
      //Tiek izveidotas reversi komplementāras sekvences IS6110 sākuma un beigu fragmentiem
      create_revcom(start_fragment_revcom, start_fragment, i);
      create_revcom(end_fragment_revcom, end_fragment, i);

      strncpy(read_start, read_sequence, i+1);
      read_start[i+1]='\0';
      strncpy(read_end, read_sequence+read_last_nuc_index-i, i+1);
      read_end[i+1]='\0';


      //printf("%s\t%d\n",read_end,i);
      //printf("%s\t%d\n",read_sequence,i);
    
      //Jāsalīdzina readu sākumi un beigas ar katru no fragmentiem
      if (strcmp(end_fragment, read_start)==0){
        strncpy(read_end_good, read_sequence+i+1,read_last_nuc_index+1);
        read_end_good[read_last_nuc_index-i]='\0';
       // OK 
       // Ieraksta pie beigām
       fprintf(fe,">end_%d\n%s\n", m,read_end_good);
       ++m;
      } 
      if (strcmp(start_fragment, read_end)==0){
        strncpy(read_start_good, read_sequence, (read_last_nuc_index-i));
        read_start_good[read_last_nuc_index-i]='\0';
       // OK
       // Ieraksta pie sākumiem
       //inplace_reverse(read_start_good);
       fprintf(fs,">start_%d\n%s\n", n,read_start_good);
       ++n;
      }
      if (strcmp(start_fragment_revcom, read_start)==0){
        strncpy(read_start_good_revcom,read_sequence+i+1,read_last_nuc_index+1);
        read_start_good_revcom[read_last_nuc_index-i]='\0'; 
        create_revcom(rev_rev_start, read_start_good_revcom, strlen(read_start_good_revcom)-1);
        //Ieraksta pie sākumiem
        //Sākumi tiek apgriezti lai IS atrastos kreisajā pusē un varētu sasortēt
        //inplace_reverse(rev_rev_start);
        fprintf(fs,">start_%d\n%s\n", n,rev_rev_start);
        ++n;
      }
      if (strcmp(end_fragment_revcom, read_end)==0){
        strncpy(read_end_good_revcom, read_sequence,(read_last_nuc_index-i));
        read_end_good_revcom[read_last_nuc_index-i]='\0';
        create_revcom(rev_rev_end, read_end_good_revcom, strlen(read_end_good_revcom)-1);
        //printf("Reads:\n%s\n\nIS end_revcom:\n%s\n\nRead end good revcom:\n%s\n\nrevrev_end:\n%s\n\n\n", read_sequence,end_fragment_revcom, read_end_good_revcom, rev_rev_end);
        //Ieraksta pie beigām
        fprintf(fe, ">end_%d\n%s\n", m, rev_rev_end);       
        ++m;
      } 
    }
  }
  //Atbrīvo funkciju malloc'd un realloc'd memory
  free(end_fragment);
  free(start_fragment_revcom);
  free(end_fragment_revcom);
  free(start_fragment);
  free(sequence);
  free(read_start_good);
  free(read_start_good_revcom);
  free(read_end_good);
  free(read_end_good_revcom);
  free(read_start);
  free(read_end);
  free(read_sequence);
  free(rev_rev_start);
  free(rev_rev_end);
  free(is_name);
  fclose(fe);
  fclose(fs);
  return 0;
}

