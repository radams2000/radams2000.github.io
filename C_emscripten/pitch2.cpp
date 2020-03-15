#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#include "signal.h"


#define FFTLEN 2048
#define FFTLEN2 1024
#define HOP 1024 // not used
#define MAXBIN 216 // about 5KHz top bin for harmonic analysis
#define MAXBLOCK 1024
#define MINBIN 40
#define FHILIMIT 1500.0
#define FLOLIMIT 25.0
//#define CONLY 1

void fft(int,int);

void hannCompute(void);
void hannMult(void);
int fftmax(void) ;
float compute_fftmagdb(void);
void Array_sort(float *, int);
void init_fft_twiddles(void);
void init_note_limits(void);
void init_acorl_tilt(void);



#ifdef CONLY
void parseWavHeader();
float get_wav16(void);
int processAudioData(float *,int, int);
FILE *wav_infile; 
FILE *infile;
FILE *outfile;

char inbuf[300];
char fname[60];
char fname_r[40];
char fname_l[40];
unsigned char ch;
int bits_per_sample,num_samples;
int num_channels,block_align;
#endif

int num_points;
int numcycles;
int numpeaks_GL = 0;

float fs_wav;
float re[FFTLEN] = {0.0};
float im[FFTLEN] = {0.0};
// float re2[FFTLEN] = {0.0};
// float im2[FFTLEN] = {0.0};
float fftpow[FFTLEN] = {0.0};

float twid_re[FFTLEN2] = {0.0};
float twid_im[FFTLEN2] = {0.0};
float acorl_tilt[FFTLEN2] = {0.0};

float fftmagdb[FFTLEN] = {0.0};
float fftpower_GL[FFTLEN] = {0.0};


float maxampl_db_GL = 0.0;
float freqs_GL[MAXBLOCK] = {0.0};



float win[FFTLEN] = {0.0};
float hann[FFTLEN] = {0.0};

float PI = 3.141592653589793;


float scale_bins2freq_GL = 1.0;
//float reliability_GL[MAXBLOCK] = {0};
int numblocks_GL;
int init = 1;
float ms_ampl_GL = 0.0;
int block_count_GL;
int num_notes_GL = 50;
float note_low_limit[50] = {0.0};
float note_nom_limit[50] = {0.0};
float note_high_limit[50] = {0.0};
int note_index_GL = 0;
float note_error_GL = 0.0;





#ifdef CONLY

int main() 
{
int i,k;
int num_return_samples = 1024;
float wavin[1000000] = {0.0};

int numblocks,block,count,tempi;
float left_in,right_in,wav_in;

outfile = fopen("freqout.txt","w");
char buff[] = "./wavs/rwa_runyon_a.wav";
if( (wav_infile = fopen(buff,"rb")) == NULL) {
	printf("cant find it\n");
	exit(0);
	}
parseWavHeader(); /** sets global num_samples, num_channels, bits_per_sample **/
printf("parsed wav file, num_samples = %d\n",num_samples);

/*********************start main routine ***************/
count = 0;
block_count_GL = 0;
hannCompute();
for(i=0;i < num_samples;i++) {

	if(num_channels == 1) { // mono
		wavin[i] = get_wav16();
	} else { // stereo
		wavin[i] = get_wav16();
	    wavin[i] = get_wav16();
	}

	

  } /** end of while loop ***/

numblocks = processAudioData(wavin,num_samples,num_return_samples);	

for(k=1;k < numblocks;k++) {
	fprintf(outfile,"%lf\n",freqs_GL[k]);

}
fclose(outfile);
}
#endif



/****************************************************/

#ifndef CONLY
extern "C" {
#endif

int processAudioData(float * wav_in, int num_samples,int num_return_samples,float fs) // returns a pointer to an array
{

int i,k,kk,count,outbuff_ptr=0;
float rtio = 0.0,ms_sum = 0.0,gain;
float max_acorl = 0.0,max_acorl2=0.0;
int argmax_acorl = 0,argmax_acorl2 =0;
float Beta,Alpha,Gamma,frac_bin,argmax_acorl_interp,raw_freq,binfund;
int num_samples_wav = num_samples-num_return_samples,maxharm,ptr,starti;
//float fs = 44.1e3; // temporary, will need a way to pass this in from the JS side
float fft_energy_1x,fft_energy_2x,fft_energy_halfx,gate,acorl_ratio,mean,freqout,tempr,pass;
int check_bin1,check_bin2,fftbin,acorlbin;
float ampl0,ampl1,ampl2,ratio_av,gate1,acorl0,acorl1,acorl_peak_pwr,acorl_nolag_pwr,kr;
float crappiness = 1.0,re_wtilt;


if(init==1) {
	init_fft_twiddles();
	init_note_limits();
	init_acorl_tilt();
	hannCompute();
	init = 0;
}



//int num_samples_wav = num_samples;
if(num_return_samples > MAXBLOCK) {
	printf("error, num return samples is < MAXBLOCK\n");
	exit(0);

}

/*********************start main routine ***************/
count = 0;
outbuff_ptr = 0;
//block_count_GL = 0;
scale_bins2freq_GL = fs_wav/(float)FFTLEN;
//printf("num_samples, num_return_samples = %d %d\n",num_samples,num_return_samples);


// for(i=num_samples_wav;i < num_samples;i++) { // zero out the return buffer
// 	wav_in[i] = 0.0;
// }

numblocks_GL = num_samples/FFTLEN; // no overlap 

for(i=0;i < num_samples;i++) {

	ms_ampl_GL = 0.998*ms_ampl_GL + 0.002*wav_in[i]*wav_in[i];

	 if(count == (FFTLEN-1)) { // non-overlapping FFT's
	 // 	for(k=0;k < 3;k++) {
		// 	printf("%2.10f  ",wav_in[k]);
		// }
	 	// load fft buffer fromm input buffer



	 	gate = 1.0;

	 	for(k=0; k < FFTLEN;k++) {re[k]=wav_in[i-FFTLEN+k+1];im[k] = 0.0;}
	 	// normalize rms level
	 	for(k=0; k < FFTLEN;k++) {ms_sum = ms_sum + re[k]*re[k];}
	 	gain = 1.0/sqrt(ms_sum/(float)FFTLEN);
	 	for(k=0; k < FFTLEN;k++) {re[k]=gain*re[k];}


		hannMult(); // in place hann window on global re
		fft(FFTLEN,1); // in-place fft on global buffer variables


		for(k=0; k < FFTLEN;k++) {
			fftpow[k] = re[k]*re[k] + im[k]*im[k];
			re[k] = fftpow[k];
			im[k] = 0.0;
		}

		fft(FFTLEN,-1);

		// zero out the early values until the autocorl goes negative
		acorl0 = re[0]; // save the 0th acorl so we can test for periodicity later
		acorl1 = re[1];
		pass = 0.0;
		for(k=0;k < FFTLEN2;k++) {
			if(re[k] < 0.0) {
				pass=1.0;
			}
			re[k] = pass*re[k];
		}
	
		// now correlt for the window tilt
		for(k=0;k < FFTLEN2;k++) {
			re[k] = re[k]*acorl_tilt[k];
		}
		// now find the peak acorl


		max_acorl = 0.0;
		argmax_acorl = 1;
		for(k=0;k < FFTLEN2;k++) { 
			if(re[k] > max_acorl) {
				max_acorl = re[k];
				argmax_acorl = k;
			}
		}

		if(argmax_acorl > 0) {
			Beta = re[argmax_acorl];
        	Alpha = re[argmax_acorl-1];
        	Gamma = re[argmax_acorl+1];
        	frac_bin = 0.5*(Alpha-Gamma)/(Alpha-2.0*Beta + Gamma);
        	argmax_acorl_interp = (float)argmax_acorl + frac_bin;
        }
        //printf("%f\n",argmax_acorl_interp);

        //check in the freq domain to see if lower-freq peaks appear
        // fftbin = (int)(((float)FFTLEN)/argmax_acorl_interp + 0.5); // fundamental should be here
        // check_bin1 = (int)(0.5*((float)FFTLEN)/argmax_acorl_interp + 0.5); // 1/2 fundamental should be here
        // check_bin2 = (int)((1.0/3.0)*((float)FFTLEN)/argmax_acorl_interp + 0.5); // 1/3 fundamental should be here
        // ampl0 = fftpow[fftbin];
        // ampl1 = fftpow[check_bin1];
        // ampl2 = fftpow[check_bin2];

        // a perfectly periodic signal will have the peak acorl = acorl[0]. If other signals are present
        // the the peak acorl wil be < acorl[0], because it never aligns 100% with any shifted version of itself
        // maybe I should modify this using the interpolation method to get the true peak instead of adding neighboring
        // bins which could possibly be negative??
        if(argmax_acorl > 0) {
        	acorl_peak_pwr = re[argmax_acorl] + re[argmax_acorl-1] + re[argmax_acorl+1];
        	acorl_nolag_pwr = acorl0 + 2.0*acorl1;
    	}

    	crappiness = acorl_nolag_pwr/acorl_peak_pwr;

        //printf("acorl_argmax, sqrt(acorl(0)), sqrt(acorl_peak_power) = %d %f %f\n",argmax_acorl,sqrt(acorl_nolag_pwr),sqrt(acorl_peak_pwr));
        // if( (fftbin > 4) && (fftbin < (FFTLEN2-4))) {
        // 	ratio_av = 0.5*(ampl0/fftpow[fftbin+3] + ampl0/fftpow[fftbin-3]);
        // 	printf("ratio_av = %f\n",ratio_av);
    
        // }
 
   		//printf("fft indexes = %d %d %d, ampl = %f %f %f\n",fftbin,check_bin1,check_bin2,ampl0,ampl1,ampl2);
        // printf("possible octave error at block %d\n",block_count_GL);
        // 	printf("comparing bin %d with bin %d\n",argmax_acorl,check_bin1);
        // 	printf("max_acorl,acorl_bin1,ratio = %f %f %f\n",max_acorl,acorl_bin1,max_acorl/acorl_bin1);
       
        //if(max_acorl/acorl_bin2 < 1.1) printf("possible triplet error at block %d\n",block_count_GL);




		//raw_freq = fs/argmax_acorl_interp;
        // octave error check, find the next-lowest-energy acorl peak
  //       if((argmax_acorl > 4) && (argmax_acorl < FFTLEN2-5)) {
  //       	// zero out the main peaks
  //       	//re[argmax_acorl] = 0;re[argmax_acorl-1] = 0;re[argmax_acorl-2] = 0;re[argmax_acorl-3] = 0;re[argmax_acorl-4] = 0;
  //       	//re[argmax_acorl+1] = 0;re[argmax_acorl+2] = 0;re[argmax_acorl+3] = 0;re[argmax_acorl+4] = 0;
		// 	starti = (int)(1.5*argmax_acorl_interp);
		// 	max_acorl2 = 0.0;
		// 	argmax_acorl2 = 1;
		// 	for(k=starti;k < FFTLEN2;k++) {
		// 		if(re[k]*acorl_tilt[k] > max_acorl2) {
		// 			max_acorl2 = re[k]*acorl_tilt[k];
		// 			argmax_acorl2 = k;
		// 		}
		// 	}
		// }
		// // if ratio is too close to 1, stop updating display
		// acorl_ratio = max_acorl/max_acorl2;
		// if(acorl_ratio < 1.1) gate = 0.0;

		//if(block_count_GL == 10) 
		// if((block_count_GL % 10)==0) {
		// 	printf("%d %f %d %f %d\n",block_count_GL,max_acorl,argmax_acorl,max_acorl2,argmax_acorl2);
		// }






        // sum the bin energies for 3 candidates, (acorl peak, acorl peak X2 and /2)
  //       if(raw_freq < 4e3) {
  //       	binfund = (float)FFTLEN/argmax_acorl_interp; // fractional bin
		// 	maxharm = (int)floor(3e3/raw_freq);
		// 	fft_energy_1x =0.0; fft_energy_2x = 0.0; fft_energy_halfx = 0.0;
		// 	for(k=1;k <= maxharm;k++) {
		// 		ptr = (int)((float)k*binfund + 0.5);
		// 		fft_energy_1x += fftpow[ptr];
		// 		ptr = (int)(2.0*(float)k*binfund + 0.5);
		// 		fft_energy_2x += fftpow[ptr];
		// 		ptr = (int)(0.5*(float)k*binfund + 0.5);
		// 		fft_energy_halfx += fftpow[ptr];

		// 	}
		// }

		// if((fft_energy_2x > fft_energy_1x) || (fft_energy_halfx > fft_energy_1x)) {
		// 	gate = 0.0;
		// 	//printf("energy wrong\n");
		// }

		// ************ This is the return to JS side!! ******************
        // if(ratio_av > 50.0) {
        // 	gate1 = 1.0;
        // } else {
        // 	gate1 = 0.0;
        // }
        gate1 = 1.0;
        if(crappiness > 1.2) gate1 = 0.0;
		freqout = gate1*fs/argmax_acorl_interp;
		if(freqout > FHILIMIT) freqout = FHILIMIT;
		//if(freqout < FLOLIMIT) freqout = FLOLIMIT;

		wav_in[outbuff_ptr] = freqout; // return pitch in same array, over-write the values that we already used 
	 	outbuff_ptr++;
	 	wav_in[outbuff_ptr] = (float)note_index_GL;
	 	outbuff_ptr++;
	 	wav_in[outbuff_ptr] = note_error_GL;



	 	// the js side will call the C program every 4K samples, so we do 2 blocks of 2K and return them in
	 	// wav_in[0] and wav_in[1]
	 	// ***********************************************************
		//printf("freq = %f\n",freqs_GL[block_count_GL]);

	 	count = 0;
	 	block_count_GL = block_count_GL + 1;
	 	outbuff_ptr++;

	} else {
		count = count + 1;
	}
	//wav_in[i] = 1.0; // debug to see if I can access this on the js side

  } /** end of for loop ***/


//for(k=0:k < 16;k++) wav_in[k] = 1;
//printf("done FFT loop, max block count =   %d\n",block_count_GL);

return block_count_GL; 

}

#ifndef CONLY
}
#endif




/**************************************************/
void hannCompute(void) {

  	int k,kk;


  	for(k = 0;k < FFTLEN;k++) {
  		win[k] = ( 0.5 * (1.0 - cos (2.0*PI*(float)k/(float)(FFTLEN-1))) );
  	}


 //  	for(k = 0;k < FFTLEN/8;k++) {
 //  		hann[k] = ( 0.5 * (1.0 - cos (2.0*PI*(float)k/(float)(FFTLEN/8-1))) );
 //  	}
 //  	for(k=0;k < FFTLEN/16;k++) {
 //  		win[k] = hann[k];
 //  	}
 //  	for(k=FFTLEN/16;k < FFTLEN-(FFTLEN/16);k++) {
	// 	win[k] = 1.0;
	// }
	// kk=0;
	// for(k=FFTLEN-(FFTLEN/16);k < FFTLEN;k++) {
	// 	win[k] = hann[FFTLEN/16+kk];
	// 	kk++;
	// }
	// for(k=0;k < FFTLEN;k++) {
	// 	printf("%f\n",win[k]);
	// }
}



/**************************************************/
void hannMult(void) {

  	int k;
  	for(k = 0;k < FFTLEN;k++) {
  		re[k] = re[k]*win[k];
  	}
}


/**************************************************/
float compute_fftmagdb(void) { // operates on global FFT registers
	int k;
	float temp;
	for(k = 0;k < MAXBIN;k++) { // only look up to 5KHz
		temp = re[k]*re[k] + im[k]*im[k];
		fftpower_GL[k] = temp;
		fftmagdb[k] = 10.0*log10(temp + 1e-12);
		if(fftmagdb[k] > maxampl_db_GL) {
			maxampl_db_GL = fftmagdb[k];
		}
		
	}
	return(maxampl_db_GL);
}
	



void init_fft_twiddles(void)
{
int k;
	
for(k=0;k < FFTLEN2;k++) {
	twid_re[k] = cos(2.0*PI * (float)(k) / (float)FFTLEN);
	twid_im[k] = sin(2.0*PI * (float)(k) / (float)FFTLEN);

}


}
/**************************************************/
void fft(int fn, int fr) // fr is 1 for forward and -1 for inverse
//int fn, ff;
/*** re and im arrays are global, so they can be huge!! **/

/* fft(re, im fn, ff) performs Fast Fourier Transform.        */
/*  re[] and im[] contain the real and imaginary part of data  */
/*   ( After transform they contain the transformed data)        */
/*  fn is the number of point. ( must be power of 2 )            */
/*  ff = -1 for forward transform.       			 */
/*     =  1 for inverse transform.       			 */

{
	float tempr, tempi, fwr, fwi, ftheta;
	//float testfwr,testfwi;
	int fj, fi, fm, mmax, tmax, istep;


	/*  do Bit-Reversal */
	fj = 1;
	for ( fi =1; fi <=fn; fi++)
	{
		if(fi < fj)
		{
			tempr = re[fj-1];
			tempi = im[fj-1];
			re[fj-1] = re[fi-1];
			im[fj-1] = im[fi-1];
			re[fi-1] = tempr;
			im[fi-1] = tempi;
		}
		fm = fn / 2.0 ;
		while ( fj > fm )
		{
			fj = fj - fm;
			fm = (fm + 1 ) / 2;
		}
		fj = fj + fm;
	}

	/* Radix two frequency decimation algorithm  */
	mmax = 1;
	tmax = FFTLEN2;
		while ( mmax < fn )
	{
		istep = 2 * mmax;
		for ( fm = 1; fm <= mmax; fm++ )
		{

			fwr = twid_re[(fm-1)*tmax];
			fwi = -((float)fr)*twid_im[(fm-1)*tmax];
			//ftheta = PI * (float)(-(fm - 1)) / (float)mmax;
			//fwr = cos(ftheta);
			//fwi = sin(ftheta);
			//printf("%d %f %f\n",(fm-1)*tmax,fwi,testfwi);
			for(fi = fm; fi <= fn; fi = fi + istep )
			{
				fj = fi + mmax;
				tempr = fwr * re[fj-1] - fwi * im[fj-1];
				tempi = fwr * im[fj-1] + fwi * re[fj-1];
				re[fj-1] = re[fi-1] - tempr;
				im[fj-1] = im[fi-1] - tempi;
				re[fi-1] = re[fi-1] + tempr;
				im[fi-1] = im[fi-1] + tempi;
			}
		}
		mmax = istep;
		tmax = tmax/2;
	}
	
}



// function to sort the array in ascending order
void Array_sort(float *array , int n)
{ 
    // declare some local variables
    int i=0 , j=0 ;
    float temp=0.0;

    for(i=0 ; i<n ; i++)
    {
        for(j=0 ; j<n-1 ; j++)
        {
            if(array[j]>array[j+1])
            {
                temp        = array[j];
                array[j]    = array[j+1];
                array[j+1]  = temp;
            }
        }
    }

    // printf("\nThe array after sorting is..\n");
    // for(i=0 ; i<n ; i++)
    // {
    //     printf("\narray_1[%d] : %d",i,array[i]);
    // }
}






/**************************************************/
#ifdef CONLY

float get_wav16()
{
int  i;
unsigned char templo,temphi;
//int fill;
int sum;
sum=0;
/**** 16-bit ************/

		templo = fgetc(wav_infile); // low byte
		if(feof(wav_infile)) printf("End 2\n");
		temphi = fgetc(wav_infile); // hi byte
		if(feof(wav_infile)) printf("End 3\n");
		//fill = 0;
		//if((temphi & 0x80) == 0x80) fill = 0xffff0000;
		//sum = ((temphi << 8) | templo) & 0x0000ffff;
		//sum = ((temphi << 8) | templo);
		temphi = (temphi + 128) & 0xffff; // invert sign bit
		sum = 256*(long)temphi + (long)templo;
		//sum = sum - 32768;
		//sum = sum | fill;
return((float)sum - 32768.0);
}
#endif



/**************************************************/
#ifdef CONLY

void parseWavHeader()
{
	/** note, format is LS byte followed by MS byte **/
int i;
int b0,b1,b2,b3,fs0,fs1,fs2,fs3,num_bytes;
for(i=1;i <= 44;i++) { 
	ch = fgetc(wav_infile); // 8 bits
	//printf("wavfile header %d %d\n",i,ch);

	if(i == 25) {fs0 = (int )ch;}
	if(i == 26) {fs1 = (int )ch;}
	if(i == 27) {fs2 = (int )ch;}
	if(i == 28) {fs3 = (int )ch;}


	if(i == 41) {b0 = (int )ch;}
	if(i == 42) {b1 = (int )ch;}
	if(i == 43) {b2 = (int )ch;}
	if(i == 44) {b3 = (int )ch;}

	if(i==23) {
		if((int)ch == 1) {
			num_channels = 1;
		}
		else {
			num_channels = 2;
		}
	}

	if(i== 35) {
		bits_per_sample = (int )ch;
	}

	}
	num_bytes = ((b3 << 24) | (b2 << 16) | (b1 << 8) | b0)/num_channels;
	if(bits_per_sample == 16) {
		num_samples = (int)(num_bytes/(2*num_channels));
	} else { // 24 bit assume
		num_samples = (int)(num_bytes/(3*num_channels));
		} 

	fs_wav = (float)(((fs3 << 24) | (fs2 << 16) | (fs1 << 8) | fs0));

	printf("fs_wav = %f\n",fs_wav);
	printf("Number of bytes in data section = %d\n",num_bytes);
	printf("Number of samples in data section = %d\n",num_samples);
	printf("num_chan= %d\n",num_channels);
	printf("bits-per-sample = %d\n",bits_per_sample);
}
#endif

void init_note_limits(void) 
	{
	int k;
	float lowest_note_freq = 440.0/8.0;// 3 octaves below A440
	for(k=0;k < num_notes_GL;k++) {
		note_nom_limit[k] = lowest_note_freq*pow(2.0,((float)k)/12.0);
		note_high_limit[k] = note_nom_limit[k]*pow(2.0,(1.0/24.0)); // 1/2 semi-tone
		note_low_limit[k] = note_nom_limit[k]*pow(2.0,(-1.0/24.0)); // 1/2 semi-tone


	} 
}

void init_acorl_tilt(void)
	{
	int k;
	double x;
	double p1,p2,p3,p4,p5,p6;
	// from polyfit of conv(hann(2048),hann(2048))
	p1 = 1.35215711132677568389e-14;
	p2 = -2.26731347185454937702e-11;
	p3 = 1.63395621017506620316e-08;
	p4 = -3.26118346124450944550e-06;
	p5 = 5.45531638598839314219e-04;
	p6 = 9.85439826782543715211e-01;

	for(k=0;k < 441;k++) { // lag 441 corresponds to a freq of 100Hz, don't go lower or too sensitive
		x = (double)(k + 1);
		acorl_tilt[k] = (float)(p6 + x*p5 + pow(x,2)*p4 + pow(x,3)*p3 + pow(x,4)*p2 + pow(x,5)*p1);
		// max value of 6.0 at k==1024
		//printf("acorltilt %d %f\n",k+1,acorl_tilt[k]);
		//if(k <= MINBIN) acorl_tilt[k] = pow((float)MINBIN,-0.25); else acorl_tilt[k] = pow( ((float)k),-0.25);
		// if(k <= 100) {
		// 	acorl_tilt[k] = pow(100.0,-0.25);
		// 	} else {
		// 		acorl_tilt[k] = pow((float)k,-0.25);
		// 	}

		//printf("acorl tilt = %f\n",acorl_tilt[k]);
	}
	for(k=442;k < FFTLEN2;k++) {
		acorl_tilt[k] =  0.0;
	}

	}





