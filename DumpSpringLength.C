
// version 4:   dump multiple layers --> points
//		dump points --> spring length, 
//		dump spring length --> spring force, according to parameters


// version 3: dump multiple layers with a list of frames 05/23/2013 >>
//	this is for the VProject file to quantify the deformation and area change	
//	also add function to compute the axial spring_length(averged over each level)
// version 2: the start node_id may begin from non-zero, 04/15/2013 >>
//		we may need to deal with multiple structures, <<
// version 1:created 04/05/2013 to dump surface node
/// STD IO
#include <iostream>
#include <sstream>
#include <fstream>
/// MY_APP
#include "CPts.h"
#include "DumpPts.h"
using namespace MY_APP;
using namespace std;
/// define area functions 
///1) at each level,we calculate area projected in x-y plane by 64-triangles
inline void 
triangle_area_fcn(const std::vector<CPts>& Level_points, const int Num_points,
		  std::vector<double>& para, int & Num_para)
{ //1) calulate the center of level, assume x=y=0;, we only need to compute z;
   // as we need to compute max, min-distance, and level cods
  // get the level and center of the area: average z coordinate
  int it_pt;
  double level_z=0.0;
  for(it_pt=0; it_pt<Num_points; it_pt++)
     level_z+=Level_points[it_pt].z;
  level_z=level_z/Num_points; // get averge 
  CPts	center(0,0,level_z);
  //CPts unitz(0,0,1.0);
  int max_id;
  int min_id;
  CPts max_pt;
  CPts min_pt;
  double max_dist=0.0;
  double min_dist=dist(Level_points[0],center);
  double area=0.0;
  //
  for (it_pt=0;it_pt<Num_points;it_pt++)
  {
    double cur_dist=dist(Level_points[it_pt],center);
    if(cur_dist>max_dist)
    { max_dist=cur_dist;
      max_id=it_pt;
      max_pt=Level_points[it_pt];
    }
    if (cur_dist<min_dist)
    {
      min_dist=cur_dist;
      min_id=it_pt;
      min_pt=Level_points[it_pt];
    }
    // compute triangle_area by cross-product
    CPts vec1=Level_points[it_pt]-center;
    CPts vec2=Level_points[((it_pt+1) % Num_points)]-center;
    double sign_area=cross_prod(vec1,vec2).z*0.5; // small triangle_area
    area+=abs(sign_area);
    
  }
  //
  //double area=PI*min_dist*max_dist;
  Num_para=6;
  para.resize(Num_para);
  para[0]=area;
  para[1]=level_z;
  para[2]=min_dist;
  para[3]=max_dist;
  para[4]=min_id;
  para[5]=max_id;
  return;
  
}


int main(int argc, char* arg[])
{// argc =2: we need to know the input file:
string InpFileInfo="mesh.info";
if (argc>1)
{ InpFileInfo=arg[1];
  cout<<"input mesh file:" << InpFileInfo << endl;
}
DumpPts my_dump(InpFileInfo);
// register new area functions
std::string default_fcn_info="area;z_level;min_dist;max_dist;min_nodeId;max_nodeId";
my_dump.RegisterComputAreaFcn(1,&triangle_area_fcn,default_fcn_info,&default_volume_fcn);
//
int Total_files=my_dump.GetTotalNum();
int Total_layers=my_dump.GetLayerNum();

string Output_File_prefix=my_dump.GetOuputPrefix();
vector< vector<double> > vec_volume;
vec_volume.resize(Total_files); // still a vector
int it_file;
for(it_file=0;it_file<Total_files; it_file++)
{
  if(!my_dump.DumpPtsBot2Top(it_file))
  { cout<< "fail to dump pts, file_no: "<<it_file<<endl;
    return -1;
  }
  vec_volume[it_file].resize(Total_layers); 
  if(!my_dump.DumpAreaBot2Top(it_file,vec_volume[it_file]))
  { cout<<"fail to dump area, file_no: "<<it_file<<endl;
    
  }
}
cout<<"++++++++++++++++++++++++++++++++++++++++++++++" << endl;
cout<<"Done: dump all the files: "<< Total_files <<endl;

// >> change 05/23/2013, consider multiple layers -> multiple volumes
string oupFileVolume=Output_File_prefix+".volume";

std::ofstream file_stream;
  file_stream.open(oupFileVolume.c_str(),std::ios::out);
  ////2.2.1) dump the file infomation
if(file_stream.is_open())
{ // 1) row dump file information
  file_stream << Total_files << "\t" << "# Number of dump timesnaps for volume" << "\n";
  // 2) dump each number per row, so we can read it easily
  for(it_file=0;it_file<Total_files;it_file++)
  {
    int it_layer=0;
    for (it_layer=0; it_layer<Total_layers; it_layer++)
    {file_stream<<vec_volume[it_file][it_layer]<<"\t";
    }
    file_stream<<"\n";
  }
  
  file_stream.close();
}
cout<<"Done: dump the volume for different files" << endl;
return 0;
 
}