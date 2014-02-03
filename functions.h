extern char * read_is_element(void);
extern void get_fragment(int start, int end,char * frag, char * seq);
extern void create_revcom(char * revcom, char * sequence, int last_index);
extern char * extract_frag(char *sequence, int is_start, int coord);
void inplace_reverse(char * str);
