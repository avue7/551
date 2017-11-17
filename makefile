tmin:
	gcc -g -Wall -W -Werror -std=c99 -o tmin tmin.c -lm

clean:
	rm -f *.o tmin
