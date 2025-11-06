# FFT implementations

I wrote this last summer because I think I wanted to do some audio visualisation in cpp but never got around to it. These are the implementations I made while I was learning about the FFT. This is also my second time ever writing cpp so I was learning alot about the language.

## Cooley turkey:

- divide and conquer recursive approach for input sizes that are power of two
- also includes the inverse using the conjugate trick

## Bluestein:

any input size, and great for prime length sequences. We use a cool identity:

>$kn = k^2 + n^2 - (k-n)^2/2$

And we are using chirp signals to transform the DFT into a convolution:

>$a[k] * exp(-i{\pi}k^2/n)$
















Putting this on my github because I like it and I want to replace one of my pinned projects.
