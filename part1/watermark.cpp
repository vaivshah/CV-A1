//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2018
//
// Ramji Chandrasekaran(ramchand)
//
//

//Link to the header file
#include <ctime>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string>
#include "SImage.h"
#include "SImageIO.h"
#include "fft.h"

using namespace std;

// function declarations
bool isPowerOfTwo(int num);
bool isImageValid(const SDoublePlane &input_imag);
SDoublePlane padImage(const SDoublePlane &input_imag);

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
  fft_real = input;
  fft_imag = SDoublePlane(input.rows(), input.cols());

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
  output_real = input_real;
  SDoublePlane output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &input_imag, const SDoublePlane &fft_real, const SDoublePlane &fft_imag) {
  SDoublePlane output_image(input_imag.rows(), input_imag.cols());
  for(int i=0; i<input_imag.rows(); i++) {
    for(int j=0; j<input_imag.cols(); j++) {
      output_image[i][j] = log(sqrt(pow(fft_real[i][j], 2) + pow(fft_imag[i][j], 2)));
    }
  }
  return output_image;
}

// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input, SDoublePlane &realPart, SDoublePlane &imaginaryPart) {
    //   SDoublePlane modifiedRealPart(input.rows(), input.cols());
    float topColAvg[realPart.cols()];
    float botColAvg[realPart.cols()];

    float topColImAvg[imaginaryPart.cols()];
    float botColImAvg[imaginaryPart.cols()];

    // filter out noise from real and imaginary buffer
    // the noise seems to be of the shape "HI" at top and bottom of the image
    // it can removed by figuring out the pixel location and setting it to the column average

    // rows 156-160, 352-356 have higher row sums - these could correspond to noise
    // calculate average of top 156 rows 
    for(int j=0; j<realPart.cols(); j++) {
        topColAvg[j] = 0;
        topColImAvg[j] = 0;
        for(int i=0; i<156; i++) {
            topColAvg[j] += realPart[i][j];
            topColImAvg[j] += imaginaryPart[i][j];
        }
        topColAvg[j] /= 156;
        topColImAvg[j] /= 156;
    }

    // calculate average of bot 154 rows 
    for(int j=0; j<realPart.cols(); j++) {
        botColAvg[j] = 0;
        botColImAvg[j] = 0;
        for(int i=357; i<realPart.rows(); i++) {
            botColAvg[j] += realPart[i][j];
            botColImAvg[j] += imaginaryPart[i][j];
        }
        botColAvg[j] /= 154;
        botColImAvg[j] /= 154;
    }

    // replace noisy rows with column averages
    for(int i=156; i<161; i++) {
        for(int j=0; j<realPart.cols(); j++) {
            realPart[i][j] = botColAvg[j];
            imaginaryPart[i][j] = topColAvg[j];
        }
    }

    for(int i=352; i<357; i++) {
        for(int j=0; j<realPart.cols(); j++) {
            realPart[i][j] = botColAvg[j];
            imaginaryPart[i][j] = botColImAvg[j];
        }
    }
    return realPart;
}

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N) {
    // declare constants to be used
    const int l = 80;
    const int r = 50;
    const float alpha = 20;
	const int midpoint = input.rows()/2;
	int digitIndex = 0;
    
    SDoublePlane realPart, imaginaryPart, output_image(input.rows(), input.cols());
	// do fft
    fft(input, realPart, imaginaryPart);
	cout<<"FFT done"<<endl;
    
    // generate random binary digits
    int randNums[l];
    // seed the random number generator
    srand(N);
    //populate the array
    for(int i=0; i<l; i++) {
        randNums[i] = rand()%2;
    }
	cout<<"generated binary digits"<<endl;

    // add watermark
	// the assignment suggest using a circle. I will rather use a square whose side
	// equals to the suggested circle. 
	// top and bottom side of the square - 2 rows and 2r columns, 20 bins in each row
	
	for(int j=midpoint-r+2, k=midpoint-r+4; j<midpoint+r, k<midpoint+r; j+=5, k+=5) {
		// place a watermark in every 5th pixel
		// top row - place watermark starting from 2nd pixel of the square
		realPart[midpoint-r][j] += (alpha * abs(realPart[midpoint-r][j]) * randNums[digitIndex++]);
		// reflect value to bottom row - symmetry
		realPart[midpoint+r][j] = realPart[midpoint-r][j];
		
		// bottom row - place watermark from 4th pixel of square
		realPart[midpoint+r][k] += (alpha * abs(realPart[midpoint+r][k]) * randNums[digitIndex++]);
		// reflect value to top row - symmetry
		realPart[midpoint-r][k] = realPart[midpoint+r][k];
	}
	cout<<"added watermark to one side"<<endl;

	// left and right side of the square - 2 columns and 2r rows, 20 bins in each column
	for(int j=midpoint-r+2, k=midpoint-r+4; j<midpoint+r, k<midpoint+r; j+=5, k+=5) {
		// place a watermark in every 5th pixel
		// left column - place watermark starting from 2nd pixel of the square
		realPart[j][midpoint-r] += (alpha * abs(realPart[j][midpoint-r]) * randNums[digitIndex++]);
		// reflect value to bottom row - symmetry
		realPart[j][midpoint+r] = realPart[j][midpoint-r];

		// right column - place watermark from 4th pixel of square
		realPart[k][midpoint+r] += (alpha * abs(realPart[k][midpoint+r]) * randNums[digitIndex++]);
		// reflect value to top row - symmetry
		realPart[k][midpoint-r] = realPart[k][midpoint+r];
	}
	cout<<"completed adding watermark"<<endl;
	
    // do ifft
    ifft(realPart, imaginaryPart, output_image);
	cout<<"IFFT complete"<<endl;
    return output_image;
}

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N) {

}

//utility function to check if a given number is power of 2
bool isPowerOfTwo(int num) {
  return (ceil(log2(num)) == floor(log2(num)));
}

// Valiate input image - checks for squareness
bool isImageValid(const SDoublePlane &input_imag) {
  if((input_imag.rows() == input_imag.cols()) && isPowerOfTwo(input_imag.rows())) {
    return true;
  }
  return false;
}

// pad a non-square image to make it square
SDoublePlane padImage(const SDoublePlane &input_imag) {
  int rows, cols;
  //ensure rows and cols are of the same length
  if(input_imag.rows() < input_imag.cols()) {
    rows = input_imag.cols();
    cols = input_imag.cols();
  } else {
    cols = input_imag.rows();
    rows = input_imag.rows();
  }
  // rows and columns might still not be power of 2
  // iteratively move to next even number until a power of 2 is found
  if (rows % 2 !=0) {
    ++rows;
    ++cols;
  }

  for(int num=2; ;num+=2) {
    if(isPowerOfTwo(rows + num)){
      rows += num;
      cols += num;
      break;
    }
  }

  cout<<"before padding: "<<input_imag.rows()<<" rows "<<input_imag.cols()<<" columns "<<endl;
  // create new image with padded cells
  SDoublePlane paddedImage(rows, cols);
  cout<<"after padding: "<<paddedImage.rows()<<" rows "<<paddedImage.cols()<<" columns "<<endl;
  // copy original image to the new image
  for(int i=0; i<input_imag.rows(); i++) {
    for(int j=0; j<input_imag.cols(); j++) {
      paddedImage[i][j] = input_imag[i][j];
    }
  }
  return paddedImage;
}


// main function
int main(int argc, char **argv)
{
  try {

    if(argc < 4)
      {
        cout << "Insufficent number of arguments; correct usage:" << endl;
        cout << "    p2 problemID inputfile outputfile" << endl;
        return -1;
      }

    string part = argv[1];
    string inputFile = argv[2];
    string outputFile = argv[3];
    cout << "In: " << inputFile <<"  Out: " << outputFile << endl;

    SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
    // input image validation
    cout<<"Validating image\n"<<"checking if image is square"<<endl;
    bool isValid = isImageValid(input_image);
    if (isValid) {
      cout<<"Image is square, proceeding to process"<<endl;
    } else {
      cout<<"Image not square, padding image"<<endl;
      input_image = padImage(input_image);
	  }

    // individual part implementation
    if(part == "1.1") {
        cout<<"Part 1.1: creating spectograms of image"<<endl;
        SDoublePlane realPart;
        SDoublePlane imaginaryPart;
        fft(input_image, realPart, imaginaryPart);
        SDoublePlane output_image = fft_magnitude(input_image, realPart, imaginaryPart);
        
        // save spectogram to output file
        // we only have greyscale image, so set R=G=B as write_png_file() takes R,G,B values
        SImageIO::write_png_file(outputFile.c_str(), output_image, output_image, output_image);
    } else if(part == "1.2") {
        cout<<"Part 1.2: filtering out noise"<<endl;
        // convert image to frequency domain and get real, imaginary planes
        SDoublePlane realPart;
        SDoublePlane imaginaryPart;
        SDoublePlane output_image(input_image.rows(), input_image.cols());
        
        // do DFT
        fft(input_image, realPart, imaginaryPart);

        // remove noise in frequency domain and do inverse fourier transform
        remove_interference(input_image, realPart, imaginaryPart);
        
        // inverse FT
        ifft(realPart, imaginaryPart, output_image);
        
        // we only have greyscale image, so set R=G=B as write_png_file() takes R,G,B values
        SImageIO::write_png_file(outputFile.c_str(), output_image, output_image, output_image);
    } else if(part == "1.3") {
		cout<<"in here"<<endl;
        if(argc < 6)
          {
            cout << "Need 6 parameters for watermark part:" << endl;
            cout << "    p2 1.3 inputfile outputfile operation N" << endl;
            return -1;
          }
        string op(argv[4]);
		const int N = atoi(argv[5]);
        if(op == "add")
          {
            // add watermark
            cout<<"adding watermark"<<endl;

            // get watermarked image
            SDoublePlane output_image = mark_image(input_image, N);

            // we only have greyscale image, so set R=G=B as write_png_file() takes R,G,B values
            SImageIO::write_png_file(outputFile.c_str(), output_image, output_image, output_image);   

          }
        else if(op == "check")
          {
            // check watermark
          }
        else
          throw string("Bad operation!");

         }
    else
      throw string("Bad part!");

  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








