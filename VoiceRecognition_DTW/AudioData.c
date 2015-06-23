#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "AudioData.h"
/*
 * input:
	 FileNameLWAV's file name
 * return :
	fileHeader:WAV's fileHeader will fill this pointer
 * author: tbc.dengwenqi@gmail.com
   all right reserved
*/
WavFileHeaderStruct* openWAVFile(char *WAVFileName)
{
	FILE* WAVFile;
	char ChunkID[4];
	WavFileHeaderStruct *fileHeader = calloc(1,sizeof(WavFileHeaderStruct));
	if (!fileHeader) return NULL;
	
	WAVFile = fopen(WAVFileName, "rb");
	if(WAVFile == NULL)
	{
//		fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		free(fileHeader);
		return NULL;	
	}
	//RIFF Chunk
	fread(&fileHeader->RIFFChunk,sizeof(fileHeader->RIFFChunk),1,WAVFile);
	if (strncmp(fileHeader->RIFFChunk.ChunkID,"RIFF",4) ||
		strncmp(fileHeader->RIFFChunk.Type,"WAVE",4))
	{
//		fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		free(fileHeader);
		fclose(WAVFile);
		return NULL;
	}
	//format chunk 
	fread(&fileHeader->FormatChunk,sizeof(fileHeader->FormatChunk),1,WAVFile);
	if (strncmp(fileHeader->FormatChunk.ChunkID, "fmt", 3))
	{
//		fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		free(fileHeader);
		fclose(WAVFile);
		return NULL;
	}
	if(fileHeader->FormatChunk.ChunkSize == 18)
	{
		//skip Additional Information
		fseek(WAVFile,2L,SEEK_CUR);
	}
	//fact chunk,optional
	fread(ChunkID,4,1,WAVFile);
	//whether comtain fact chunk
	if(!strncmp(ChunkID,"fact",4))
	{
		
		fread(&fileHeader->FactChunk.ChunkSize,sizeof(fileHeader->FactChunk.ChunkSize),1,WAVFile);
		fread(&fileHeader->FactChunk.Data,sizeof(fileHeader->FactChunk.Data),1,WAVFile);
		fread(fileHeader->DataChunk.ChunkID,4,1,WAVFile);	
		ChunkID[0] = 'n';
	}
	else
	{
		fileHeader->FactChunk.ChunkSize = 4;
		fileHeader->FactChunk.Data = 0;
		strncpy(fileHeader->DataChunk.ChunkID, ChunkID, 4); 
	}
	if(strncmp(fileHeader->DataChunk.ChunkID,"data",4))
	{
//		fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		free(fileHeader);
		fclose(WAVFile);
		return NULL;
	}
	fread(&fileHeader->DataChunk.ChunkSize,sizeof(fileHeader->DataChunk.ChunkSize),1,WAVFile);	
	fileHeader->DataChunk.data = (uint8_t *)calloc(fileHeader->DataChunk.ChunkSize,sizeof(uint8_t));
	fread(fileHeader->DataChunk.data,1,fileHeader->DataChunk.ChunkSize,WAVFile);
	
	fclose(WAVFile);
	return fileHeader;
}
int closeWAVFIle(WavFileHeaderStruct * fileHeader)
{
	if (!fileHeader)
	{
	//	fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		return -1;
	}
	if (!fileHeader->DataChunk.data)
		free(fileHeader->DataChunk.data);
	free(fileHeader);
	return 0;
}
/*
 *  
 */
float *readAllWAVData(WavFileHeaderStruct *fileHeader,int *n)
{
	int16_t *dataInt16;
	int16_t max;
	int i;
	float *data = NULL;

	if (!fileHeader || !fileHeader->DataChunk.data)
	{
	//	fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		return NULL;
	}

	*n =  fileHeader->DataChunk.ChunkSize / (fileHeader->FormatChunk.BitsPerSample / 8);
	if (!n)
	{
	//	fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		return NULL;
	}
	data = calloc(*n, sizeof(float));
	if (!data) return NULL;

	if (fileHeader->FormatChunk.NumChannels == 2)
	{
		//双声道，暂不处理
	//	fprintf(stderr, "error：%s,%s\n", __FILE__, __LINE__);
		return NULL;
	}
	if (fileHeader->FormatChunk.BitsPerSample == 8)
	{
		max = fileHeader->DataChunk.data[0];
		for (i = 1; i < *n; ++i)
			if (fileHeader->DataChunk.data[i] > max)
				max = fileHeader->DataChunk.data[i];
		if (!max) max = 1;
		for (i = 0; i < *n; ++i)
			data[i] = fileHeader->DataChunk.data[i] / (float)max;
	}
	else if (fileHeader->FormatChunk.BitsPerSample == 16)
	{
		dataInt16 = (int16_t*)fileHeader->DataChunk.data;
		max = dataInt16[0];
		for (i = 1; i < *n; ++i)
			if(dataInt16[i] > max)
				max = dataInt16[i];
		if (!max) max = 1;
		for (i = 0; i < *n; ++i)
			data[i] = dataInt16[i] / (float)max;
	}
	return data;
}
