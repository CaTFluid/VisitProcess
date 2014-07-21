// test Discard_comments unitlity
#include <DiscardComments.h>
#include <string>
#include <iostream>
#include <fstream>
using namespace MY_CPP;
int main()
{// follow the TimeDependentSpring
  std::string line_string;
  std::ifstream file_stream;
  std::string fileName="test.log";
  // first write, or output
  std::ofstream file_out;
  file_out.open(fileName.c_str(), std::ios::out);
  if(file_out.is_open())
    {
      file_out << "test 123 //123 test";
      file_out.close();
    }
  else // print error
    std::cout<< "can not write test.log"<<std::endl;
  // second dump to test
    

  file_stream.open(fileName.c_str(), std::ios::in);
  if (file_stream.is_open())
    {
      if(!std::getline(file_stream,line_string))
	// error
	std::cout<< "error for read"<<std::endl;
      else
	{// line_string gets data
	  std::cout<<"before discard comments"<<
	    line_string<< std::endl;
	  line_string=discard_comments(line_string);
	  std::cout<<"after discard comments"<< 
	    line_string<<std::endl;
	}
      file_stream.close();
    }
}
