#include <iostream>
#include <cmath>
#include <smithlin.h>
#include <fstream>

using namespace std;

const int N = 100000; //number of particles
const int C = 100; //number of collisions per output
const int D = 10000; //number of outputs
const double pi = 3.14159265; //pi, obvs
double rms(Vector v[N]); //calculates root-mean-square speed of particle system
void write(int, Vector v[N]); //outputs particle speeds to text file
bool collide(int, int, Vector v[N]); //forces two particles to undergo a collision
int rand_part(); //custom random integer generator

int main()
{
	srand(44547998);
	ofstream uni;
	cout << "System: " << N << " particles; " << C << " collisions per output" << endl; 
	Vector v[N]; //particle velocities
	for(int i = 0; i < N; i++)
	{ //initialising particle speeds to follow uniform distibution and for direction to be randomised
		double s = rand() % 2000;
		v[i] = Vector(3);
		double theta = (double)(rand() % (int)(2000*pi))/1000.0;
		double phi = (double)(rand() % (int)(2000*pi))/1000.0;
		v[i].a[0] = s*cos(theta)*cos(phi);
		v[i].a[1] = s*sin(theta)*cos(phi);
		v[i].a[2] = s*sin(phi);
	}
	cout << "Initial v(rms) = " << rms(v) << "; a = " << rms(v)/sqrt(3) << endl; //system properties
	for(int d = 0; d < D; d++)
	{
		write(d, v); //output to text file (ensure folder 'output' exists
		for(int c = 0; c < C; c++) 
		{
			int p1 = rand_part(); //selecting particles
			int p2 = rand_part();
			while(!collide(p1, p2, v)) //colliding, whilst neglecting trivial collisions (i.e if p1 = p2 or collision impossible)
			{
				p1 = rand_part();
				p2 = rand_part();
			}
		}
	}
	cout << "Final v(rms) = " << rms(v) << "; a = " << rms(v)/sqrt(3) << endl; //systems properties; should be idetical to initial properties
	for(int i = 0; i < N; i++) v[i].destroy(); //deallocating memory
	return 0;
}
void write(int d, Vector v[N])
{
	ofstream output;
	char filename[40];
	sprintf(filename, "output/d%d.txt", d); //filename named after output number
	output.open(filename);
	for(int i = 0; i < N; i++) output << v[i].modulus() << endl; //writing
	output.close();
	return;
}
bool collide(int p1, int p2, Vector v[N])
{
	bool isCollide;
	Vector r = Vector(3);
	double theta = (double)(rand() % (int)(2000*pi))/1000.0; //randomising collision orientation
	double phi = (double)(rand() % (int)(2000*pi))/1000.0;
	r.a[0] = cos(theta)*cos(phi);
	r.a[1] = sin(theta)*cos(phi);
	r.a[2] = sin(phi);
	Vector Du = v[p1] - v[p2];
	double k = r*Du;
	Vector Im = k*r; 
	if(k <= 0) isCollide = false; //checking for triviality
	else
	{
		isCollide = true;
		v[p1] = v[p1] - Im;
		v[p2] = v[p2] + Im;
	}
	r.destroy(); //deallocating memory
	Du.destroy();
	Im.destroy();
	return isCollide;
}
double rms(Vector v[N])
{
	double s = 0;
	for(int i = 0; i < N; i++) s += v[i]*v[i];
	s /= N;
	s = sqrt(s);
	return s;
}
int rand_part()
{
	int rnd;
	int div = (int)ceil((double)N/RAND_MAX);
	retry:
	int a = rand() % div;
	if(a == div - 1)
	{
		int b = N % RAND_MAX;
		rnd = rand() % RAND_MAX;
		if(rnd >= b) goto retry;
	}else
	{
		rnd = rand() % RAND_MAX;
	}
	rnd += a*RAND_MAX;
	return rnd;
}
