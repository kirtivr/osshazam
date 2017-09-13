I implemented Shazam based on the creators' original paper
http://www.ee.columbia.edu/~dpwe/papers/Wang03-shazam.pdf.
The results have been encouraging.

How to test:
* Get the audio for the songs you want to compare.

On  Audacity:
1. Convert stereo mp3 (2 channels – meant for the right and left
speakers respectively – convey the impression of sound coming
from 2 directions) to mono(single channel) .wav format.
2. Clip the songs at the point where you think you have enough information.
3. Sampled at 16000 samples per second.
4. Export samples from audacity to a text file.
* Note that I tested with 848000 samples – in the range of 2^20.

In the file “songcompare.cpp” line 469, 470 and 471, replace the songs with the samples you want to compare. 

* You need Cilk runtime with clang to compile and run.
