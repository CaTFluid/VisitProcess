// test pts.h and ifstream
#include "CPts.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;
using namespace MY_APP;
int main()
{
string FileName="4file.log";
int num_pts=12;
ifstream file_stream;
file_stream.open(FileName.c_str(), std::ios::in);
vector<CPts> int_pts(num_pts);
int it_pt;
for (it_pt=0;it_pt<num_pts;it_pt++)
{cout << int_pts[it_pt] << endl;
}
if(file_stream.is_open())
{
 for (it_pt=0; it_pt<num_pts; it_pt ++)
	{file_stream >> int_pts[it_pt];
	}

 file_stream.close();
 
}
cout<< "after read:"<<endl;
for (it_pt=0;it_pt<num_pts;it_pt++)
{cout << int_pts[it_pt] << endl;
}
}
 
