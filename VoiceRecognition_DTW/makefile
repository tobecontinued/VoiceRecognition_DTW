cc = gcc
CFLAGS = -g
objs = AudioData.o VoiceRecognition.o MFCC.o main.o
VoiceRecognition_DTW:$(objs) 
	$(cc) $(CFLAGS) -o $@ $(objs) -lm
.c.o:$(cc) $(CFLAGS) $<

clean:
	rm -f *.o VoiceRecognition_DTW 
