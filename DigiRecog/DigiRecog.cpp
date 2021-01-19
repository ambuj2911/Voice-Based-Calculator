// DigiRecog.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "vector"
#include  "stdio.h"
#include <sstream>
#include "iostream"
#include <fstream>
#include <cmath>
#define ld long double
using namespace std;

vector<ld>pi;
vector< vector<ld> > codebook, a, b, alpha, beta, new_a, new_b;
int n = 5, m = 32;
vector<vector<vector<ld>>> si, gamma;


void codebook_gen()
{
	ifstream myfile;
	myfile.open("codebook.txt");
	ld data;
	int count = 0;
	vector<ld> temp;
	while(myfile >> data)
	{
		if(count % 12 == 0 && count != 0)
		{
			codebook.push_back(temp);
			temp.clear();
		}
		temp.push_back(data);
		count++;
	}
	codebook.push_back(temp);
}


vector<ld> Ri_Calculate(vector <ld> s, int p)
{
	 int k, i;
	 vector <ld> r;
	 for(k=0;k <= p; k++)
	 {
		 ld temp = 1;
		 ld sum = 0;
		 for(i = 0; i <= 319 - k ; i++)
		 {
		   temp = s[i] * s[i+k];
		   sum = sum + temp;
		 }
		 r.push_back(sum);
	 }
	 return r;
}



vector<ld> Ai_Calculate(vector <ld> r, int p)
{
	vector<ld> ai;
	ai.push_back(1);
	ai.push_back(0);
	int i , j, t = 1;
	ld k = (r[1] / r[0])* -1;
	ai[1] = k;
	ld a = r[0] * (1 - (k * k));
	vector <ld> temp_a;
	
	for(i = 2;i <= p;i++)
	{
		ld s = 0;
		for(j = 1; j <= i-1; j++)
		{
			s = s + (r[j] * ai[i-j]);
		}
	
		s = s + r[i];
		k = (s/a) * -1;

		temp_a.clear();
		temp_a.push_back(0);

		for(j = 1; j <= i-1; j++)
		{
			temp_a.push_back( ai[j] + k * ai[i-j]);
		}

		for(j = 1; j<= i-1;j++)
		{
			ai[j] = temp_a[j];
		}
		
		ai.push_back(k);
		a = a * (1 - (k*k));

	}

   

	for(i = 0; i< ai.size() ; i++)
	{
		ai[i] = ai[i] * -1;
	}
	
	
	return ai;

}


vector<ld> raised_sine(vector <ld> c)
{
	ld win;
	for(int i = 0; i < 12; i++)
	{
		win = 6 * sin((3.14 * (i + 1)) / 12);
		win++;
		c[i] = c[i] * win;
	}

	return c;
}


vector<ld> Ci_Calculate(vector <ld> a, int p)
{
	int i,k;
	ld temp;
	vector<ld> c;
	c.push_back(0);
	for(i = 1; i<=p; i++)
	{
		temp = a[i];
		for(k = 1; k <= i- 1; k++)
		{
			temp += (c[k]* a[i-k]*k) / i ;
		}
		c.push_back(temp);
	}
	c.erase(c.begin());
	return raised_sine(c);
}


ld tokuhara(vector<ld> ci_train, vector<ld> ci_test)
{
	ld difference,tokuhara_distance, final_distance = 0;
	ld weight_T[]={1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};
	
	tokuhara_distance=0;
	for(int j=0;j<12;j++)
	{
		difference=ci_train[j]-ci_test[j];
		tokuhara_distance+=(difference*difference*weight_T[j]);
	}
	final_distance+=tokuhara_distance;
	
	return (final_distance/5);
}



int obs_seq(vector<ld> frame)
{
	vector<ld> a,r,c;
	r = Ri_Calculate(frame, 12);
	a = Ai_Calculate(r, 12);
	c = Ci_Calculate(a,12);
	ld dist, min;
	int mark = 1;
	min = tokuhara(c, codebook[0]);
	for(int i = 1; i < codebook.size(); i++)
	{
		dist = tokuhara(c, codebook[i]);
		if(dist < min)
		{
			min = dist;
			mark = i + 1;
		}
		
	}
	return mark;
}

void print_2(vector< vector<ld> > temp)
{
	for(int i = 0; i < temp.size(); i++)
	{
		for(int j = 0; j < temp[i].size(); j++)
			cout << temp[i][j] << " ";
		cout << endl;
	}

}



void init()
{
	ifstream myfile, bfile, pifile;
	myfile.open("init_val\\A_MATRIX.txt");
	ld data;
	int count = 0;
	vector<ld> temp;
	while(myfile >> data)
	{
		if(count % 5 == 0 && count != 0)
		{
			a.push_back(temp);
			temp.clear();
		}
		temp.push_back(data);
		count++;
	}
	a.push_back(temp);
	temp.clear();
	bfile.open("init_val\\B_MATRIX.txt");
	count = 0;
	while(bfile >> data)
	{
		if(count % 32 == 0 && count != 0)
		{
			b.push_back(temp);
			temp.clear();
		}
		temp.push_back(data);
		count++;
	}
	b.push_back(temp);

	pifile.open("init_val\\PI_MATRIX.txt");
	count = 0;
	while(pifile >> data)
	{
		pi.push_back(data);
	}

	/*for(int i = 0; i < b.size(); i++)
	{
		for(int j = 0; j < b[i].size(); j++)
			cout << b[i][j] << " ";
		cout << endl;
	}*/

}


void forward(vector<int> obs, int index)
{
	if(index == obs.size())
		return;
	if(index == 0)
	{
		vector<ld> temp;
		ld prod;
		for(int i = 0; i < 5; i++)
		{
			prod = pi[0] * b[i][obs[0]];
			temp.push_back(prod);
		}
		alpha.push_back(temp);
		forward(obs, index + 1);
	}
	else
	{
		vector<ld> temp;
		for(int i = 0; i < n; i++)
		{
			ld sum = 0;
			ld prod = 1;
			for(int j = 0; j < alpha[index - 1].size(); j++)
			{
				prod = alpha[index - 1][j] * a[j][i];
				sum = sum + prod;
			}
			sum = sum * b[i][obs[index]];
			temp.push_back(sum);
		}
		alpha.push_back(temp);
		forward(obs, index + 1);
	}
	return;
}




void backward(vector<int> obs, int index)
{
	if(index == -1)
		return;
	else
	{
		for(int i = 0; i < n; i++)
		{
			ld sum = 0, prod = 1;
			for(int j = 0;j < n; j++)
			{
				prod = a[i][j] * b[j][obs[index + 1]] * beta[index + 1][j];
				sum += prod;
			}
			beta[index][i] = sum;
		}
		backward(obs, index - 1);
	}


}



void hmm_matrix(vector<int> obs)
{
	vector<ld> prev(5);
	forward(obs, 0);
	vector<ld> temp(5, 0);
	for(int i = 0; i < obs.size() - 1; i++)
		beta.push_back(temp);
	temp[0] = 1;
	temp[1] = 1;
	temp[2] = 1;
	temp[3] = 1;
	temp[4] = 1;
	beta.push_back(temp);
	
	backward(obs, obs.size() - 2);
	
	print_2(beta);
	// si matrix
	// gamma matrix

}






vector <ld> find_utterance()
{
	ifstream myfile;
	myfile.open("example.txt");
	char str[60];
	vector<long> marker;
	vector<ld> ste_e;
	vector<ld> ste_frame;
	FILE *fp = fopen("example.txt", "r");
	int k = -1;
	ld x;
	ld ste = 0;
	long int line = 0;
	vector<ld> val;
	while( fgets (str, 60, fp)!=NULL ) { // Iterating over lines of files.    
	  myfile >> x;
	  val.push_back(x);
	  if(line % 80 == 0 && line != 0)
	  {
		  ste_e.push_back(ste);
		  ste = 0;
	  }
	  ste += (x * x);
	  line++;
	}
	if(line % 80 != 0)
	{
		ste_e.push_back(ste);
		ste = 0;
	}
	ld max = -1;
	long long h;
	for(int i = 0; i < ste_e.size(); i++)
	{
		ld ste = 0;
		int j;
		for(j = 0; j < 4; j++)
		{
			if(i + j < ste_e.size())
			{
				ste += ste_e[i + j];
			}
		}
		ste /= (80 * j);
		ste_frame.push_back(ste);
		if(ste > max)
		{
			max = ste;
			h = i;
		}
	}
	// ste_frame has ste for all frames 
	// from here just for example frame taken from 0 line
	// do word starting calc here.
	vector< vector<ld> > temp;
	
	for(int i = 0; i < 85; i++)
	{
		vector<ld>result(320);
		vector<ld>::iterator ptr; 
		ptr = val.begin() + (i * 80);
		copy(ptr, ptr + 320, result.begin()); 
		ld ham_weight;
		line = 0;
		for(int i = 0; i<= 319; i++)
		{
			ham_weight = 0.54 - 0.46 * cos(2 * 3.14 * i/ 319);
			result[i] *= ham_weight;
	   }
		temp.push_back(result);
	}

	vector<int> obs;
	for(int i = 0; i < 85; i++)
	{
		int mark = obs_seq(temp[i]);
		obs.push_back(mark - 1);
	}

	// obs has observation seq start with hmm;
	hmm_matrix(obs);
	return ste_frame;


}




ld normalise(FILE * fp, ld dc)
{
	char str[60];
	ld max = 0, min = 0;    
	if(fp == NULL) {
      perror("Error opening file");
      return(-1);
      }
	
	while( fgets (str, 60, fp)!=NULL ) { // Iterating over lines of files.    
	  stringstream num(str);// Lines of file are readed as strings so a converter is used to convert stringg to integer.
      ld x = 0;
	  num >> x;   
	  x = x - dc;
	
	  if(max < x)
	  {
		max = x;
	  }
	  if (min > x)
	  {
		min = x;
	  }
	   
	}

	if (max > min * -1)
	{
	  return(10000 / max);
	}
	else
	{
	  return(10000 / (min*-1));
	}

}



void normalise_file(FILE *fp, ld norm_fact, ld dc)
{
	char str[60];
	ofstream myfile;
	long long int line = 0;
	if(fp == NULL) {
      perror("Error opening file");
      
      }
	myfile.open ("example.txt");
	if(myfile == NULL) {
      perror("Error opening file example");
      
      }
	while( fgets (str, 60, fp)!=NULL ) { // Iterating over lines of files.
    
	  stringstream num(str);// Lines of file are readed as strings so a converter is used to convert stringg to integer.
      ld x = 0;
	  num >> x;    
	  myfile << (x - dc) * norm_fact;
	  myfile<<"\n";
	  line++;
	}
	
}



ld DC_Shift(FILE *fp)
{
  
  ld dc = 0;
  char str[60];
  int line = 0, flag = 0, k = 1;
  
  while( fgets (str, 60, fp)!=NULL && k == 1 ) { // Iterating over lines of files.
    
	  stringstream num(str);// Lines of file are readed as strings so a converter is used to convert stringg to integer.
      ld x = 0;
	  num >> x;    
	  if (flag == 3)
	  {
		dc = dc / (3 * 320);
		k = 0;
	  }
	  else
	  {
		 if (line % 320 == 0 && line != 0)
		 {
			dc = dc + x;
			flag++;
			line++;
			
		 }
		 else
		 {
		   dc = dc + x;
		   line++;
		 }
	  }
	}
  
  return dc;

}



int _tmain(int argc, _TCHAR* argv[])
{
	FILE* fp;
	fp = fopen("1.txt", "r");
	ld dc_corr = 0, norm_fact;
	dc_corr = DC_Shift(fp); // DC_shift calculations.
	norm_fact = normalise(fp, dc_corr);
	fp = fopen("1.txt", "r");
	char str[60];
	normalise_file(fp, norm_fact, dc_corr);
	vector<ld> marker;
	codebook_gen();
	init();
	marker = find_utterance();
	return 0;
}

