I discovered a tmin of about 10,500,000. I decided to use a long double so there would be less error in the calculations and to make sure the calculations could actually fit to the required number of sig figs.
To find tmin I ran my code that iterates through all the possible n's starting at the given t and gives the result for each n. I ran it many times and let it print the error and result of each iteration to know how close I was.
As I increased the starting t, I got closer and at 10,500,000 for n the absolute relative true error was changing between approximately 4.9703968080889452593e-15 and 5.0129342916510645811e-15 depending on whether the n was even or odd.
Since 4.9703968080889452593e-15 is smaller than the 0.5e-14 stopping criteria, then I determined that this n of 10,500,000 must be our tmin (or close to it).

It did not take that long to find the tmin by incrementing the t by a large amount and re-running the program to see how much closer the error was to our stopping criteria.
However, if I had just started the program with a t of 1 and incremented my n by 1 each iteration, my program would probably take days to complete. Doing the trapazoid rule with an n of 10,500,000 takes over a second on my computer.
To find a tmin faster with the computer I did also modify how much the n was incremented each iteration. I changed it to only increment n by 1 though once I was getting closer to the tmin (so the program in its current state only increments n by one).
Probably, letting the program dynamically choose how much to increment n by would be a good approach, but I didn't try it.
