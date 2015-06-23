#ifndef _GET_AUDIO_DATA_H_
#define _GET_AUDIO_DATA_H_
#include <stdint.h>
///////////WAV file header////////// 
//ref:http://baike.baidu.com/view/8033.htm?fromId=14471
typedef struct
{
	char ChunkID[4]; // "RIFF"; The "RIFF" the mainchunk;
	uint32_t ChunkSize; // FileSize – 8; The size following this data
	char Type[4]; // "WAVE"; The "WAVE" format consists of two subchunks: "fmt " and "data"
}RIFFRStruct; 
typedef struct
{
	char ChunkID[4]; // "fmt "
	uint32_t ChunkSize; // 16 for PCM. This is the size of the rest of the subchunk which follows this data.
	uint16_t AudioFormat; // 1 for PCM. Linear quantization
	uint16_t NumChannels; // 1->Mono, 2->stereo, etc..
	uint32_t SampleRate; // 8000, 11025, 16000, 44100, 48000, etc..
	uint32_t ByteRate; // = SampleRate * NumChannels * BitsPerSample/8
	uint16_t BlockAlign; // = NumChannels * BitsPerSample / 8
	uint16_t BitsPerSample; // 8->8bits, 16->16bits, etc..每个样本的位数
	//unsigned short AdditionalInformation; //ChunkSize = 18 时才有
}FormatStruct;
typedef struct
{
	//optional	
	char ChunkID[4]; // "fact "
	uint32_t ChunkSize; // = 4
	uint32_t Data;
}FactStruct;
typedef struct
{
	char ChunkID[4]; // "data "
	uint32_t ChunkSize;//数据的中总共的byte数
	uint8_t *data;
}WAVDataStruct;
typedef struct 
{
	RIFFRStruct  RIFFChunk;
	FormatStruct FormatChunk;
	FactStruct FactChunk;//optional for wav file 
	WAVDataStruct DataChunk;
}WavFileHeaderStruct;	
WavFileHeaderStruct* openWAVFile(char *WAVFileName);
int closeWAVFIle(WavFileHeaderStruct * fileHeader);
float *readAllWAVData(WavFileHeaderStruct *fileHeader, int *n);
#endif
