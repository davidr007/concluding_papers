#include <stdio.h>
#include <stdlib.h>







int main(int argc, char* argv[])
{
  int tmp;
  float realtmp;
  int nread = 1;
  FILE* f = fopen(argv[1], "rb");
  int count = 0;
  int op_mode;

  
  op_mode = atoi(argv[2]); /* 1=>3int_3rea: 2=>int: 3=>real */
  
  while(1)
    {
      realtmp = 0.0;
      tmp = 0;
		
      //nread = fread(&realtmp, 4, 1, f);
      if(nread == 0)
	{
	  break;
	}

      if(op_mode == 1)  /* 3 integers then 3 reals */
	{
	  if(count % 6 < 3)
	    {
	      nread = fread(&tmp, 4 /* # bytes to read */, 1, f);
	      printf("%d: %d\n", count, tmp);
	    }
	  else
	    {
	      nread = fread(&realtmp, 4, 1, f); 
	      printf("%d: %f\n", count, realtmp); 
	    }
	}
      
      if(op_mode == 2)  /* all integers */
	{
	  nread = fread(&tmp, 4 /* # bytes to read */, 1, f);
	  printf("%d: %d\n", count, tmp);
	}

      if(op_mode == 3)
	{
	  nread = fread(&realtmp, 4, 1, f); 
	  printf("%d: %e\n", count, realtmp); 
	}

      count++;
    }
  fclose(f);	

  return 0;

}
