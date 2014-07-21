// modified: 12/19/2012:  add mpi_sum reduction to evaluate actually updated springs;
//  
// File: TimeDependentSpring, concrete class to update spring constants
#include "TimeDependentSpring.h"

///////////////////////// INCLUDES /////////////////////////

// IBAMR INCLUDES
//#include <ibamr/IBBeamForceSpec.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
//#include <ibtk/LNodeSet.h>
#include <ibtk/LNodeSetData.h>

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>
// STD C++ INCLUDES


///////////////////////// NAMESPACE ////////////////////////

namespace IBAMR
{
namespace
{

inline std::string
discard_comments(
    const std::string& input_string)
	{
    // Create a copy of the input string, but without any text following a '!',
    // '#', or '%' character.
    std::string output_string = input_string;
    std::istringstream string_stream;

    // Discard any text following a '!' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '!');
    string_stream.clear();

    // Discard any text following a '#' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '#');
    string_stream.clear();

    // Discard any text following a '%' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '%');
    string_stream.clear();
    return output_string;
	}// discard_comments
}
  /////////////////////// PUBLIC ///////////////////////////

TimeDependentSpring::TimeDependentSpring(
	Pointer<PatchHierarchy<NDIM> > hierarchy,
        LDataManager* const lag_manager,
        std::string SpringInfoFile)
    : d_hierarchy(hierarchy),
      d_lag_manager(lag_manager),
      d_SpringInfoFile(SpringInfoFile),
      d_no_mod_springs(-1),
      d_check_timeSteps(0) // by default: no check; 
	     
{
    std::string line_string;
    std::ifstream file_stream;
    std::string object_name="[update_spring]";
    file_stream.open(d_SpringInfoFile.c_str(), std::ios::in);
    if (file_stream.is_open())
    {
        pout << "++++++++++++++  See File: "<< d_SpringInfoFile << "   +++++++++++++++ \n";
        pout << d_SpringInfoFile << ":  "
             << "Processing time dependent spring information" <<d_SpringInfoFile << std::endl;


        // The first entry in the file is the number of mod_springs.
        if (!std::getline(file_stream, line_string))
        {
            TBOX_ERROR(object_name << ":\n  Premature end to input file encountered before line 1 of file " << d_SpringInfoFile << std::endl);
        }
        else
        {
            line_string = discard_comments(line_string);
            std::istringstream line_stream(line_string);
            if (!(line_stream >>d_no_mod_springs))// automatic convert type
            {
                TBOX_ERROR(object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << d_SpringInfoFile << std::endl);
            }
            pout<< d_SpringInfoFile << ":  "<<"Number of modified springs = "<< d_no_mod_springs<< "\n";
        }

        // Second entry is the number of spring info--> means how many continuous springs in (min, max) nodal range

        if (!std::getline(file_stream, line_string)) {
            TBOX_ERROR(object_name << ":\n  Premature end to input file encountered before line " << 2 << " of file " << d_SpringInfoFile << std::endl);// some announcement
        }
        else
        { line_string=discard_comments(line_string);
            std::istringstream line_stream(line_string);
            if (!(line_stream >>d_no_spring_info))// automatic convert type
            {
                TBOX_ERROR(object_name << ":\n  Invalid entry in input file encountered on line 2 of file " << d_SpringInfoFile << std::endl);
            }
        }
        pout<< d_SpringInfoFile << ":  "<<"Number of spring information = "<<d_no_spring_info<< "\n";
        // Each successive line provides the min and max and spring_fnx
        d_spring_info.resize(d_no_spring_info); // resize the vector
        d_wave_info.resize(d_no_spring_info);
        d_bound_check.resize(d_no_spring_info);// bounds of rest_length, added 03/23/2013
        for (int k = 0; k < d_no_spring_info; ++k)
        {
            if (!std::getline(file_stream, line_string))
            {
                TBOX_ERROR(object_name << ":\n  Premature end to input file encountered before line " << 3*k+2 << " of file " << d_SpringInfoFile << std::endl);
            }
            else // read data into d_spring_info;
            {
                line_string = discard_comments(line_string);
                std::istringstream line_stream(line_string);
                line_stream>>d_spring_info[k].idx_min_node;
                line_stream>>d_spring_info[k].idx_max_node;
                line_stream>>d_spring_info[k].idx_spring_fn;
                pout << "++++++++++++++				 +++++++++++++++\n";
                pout<<d_SpringInfoFile << ":  "<<"Spring information of list: " << k<< "  >> \n"
                   <<"      ++++ Global index of start node = "<< d_spring_info[k].idx_min_node << "\n"
                  <<"      ++++ Global index of end   node = "<< d_spring_info[k].idx_max_node << "\n"
                 <<"      ++++ Index of spring function   = "<< d_spring_info[k].idx_spring_fn<< "\n"
                <<d_SpringInfoFile << ":  "<<"Spring information of list: " << k<< "  << \n";
                pout << "++++++++++++++ 		   	+++++++++++++++\n"; }
            if (!std::getline(file_stream, line_string))
            {
                TBOX_ERROR(object_name << ":\n  Premature end to input file encountered before line " << 3*k+3 << " of file " << d_SpringInfoFile << std::endl);
            }
            else // read data into d_wave_info;
            {
                line_string = discard_comments(line_string);
                std::istringstream line_stream(line_string);
		// first line: index of wave algorithm ,  number of wave parameters
                line_stream>>d_wave_info[k].idx_wave_algorithm;
                line_stream>>d_wave_info[k].num_wave_parameter;
                d_wave_info[k].vec_wave_parameter.resize(d_wave_info[k].num_wave_parameter);
		
                std::getline(file_stream, line_string);
		pout << "++++++++++++++				 +++++++++++++++\n";
		pout << "	+++++ wave idx = "<< d_wave_info[k].idx_wave_algorithm <<" +++++   \n";
		pout << "	+++++ wave parameters +++++   \n";
		pout << "	+++++" << line_string << "\n";
                line_string = discard_comments(line_string);
                std::istringstream line_stream2(line_string);
		// next line: wave parameters, including tube_length, and axial_mesh_size
                for (int kpara=0; kpara<d_wave_info[k].num_wave_parameter;++kpara)
                {line_stream2>>d_wave_info[k].vec_wave_parameter[kpara];

                }
                //>> to check the wave parameters added 01/10/2013
		pout << "++++++++++++++                          +++++++++++++++\n";
		pout << "	+++++ check vec_wave +++++	\n";
		for (int kpara=0; kpara<d_wave_info[k].num_wave_parameter;++kpara)
                {pout << d_wave_info[k].vec_wave_parameter[kpara]<<"; ";

                }
		pout << "\n";
		//<< 01/10/2013
		//>> to add upper bound, lower bound for check the spring rest length: bug during regrid 03/23/2013 
		std::getline(file_stream, line_string);
		pout << "++++++++++++++				 +++++++++++++++\n";
		pout << "	+++++ lower/upper bounds of springs' rest_length +++++   \n";
		pout << "	+++++" << line_string << "\n";
		line_string = discard_comments(line_string);
                std::istringstream  line_stream1(line_string);
                line_stream1>>d_bound_check[k].lower_check;
                line_stream1>>d_bound_check[k].upper_check;
		if (d_bound_check[k].lower_check > d_bound_check[k].upper_check)
		{ double temp=d_bound_check[k].lower_check;
		  d_bound_check[k].lower_check=d_bound_check[k].upper_check;
		  d_bound_check[k].upper_check= temp;
		}
		pout << "	+++++ check lower/upper bounds +++++	\n";
		pout << d_bound_check[k].lower_check << "; " << d_bound_check[k].upper_check << "\n";
		//<< to add check bounds, 03/23/2013
            }
        }
        // Then we read information for the mesh node index_geom. structure
        // this should be one line: 3 integers.
        if (!std::getline(file_stream, line_string))
        {
            TBOX_ERROR(object_name << ":\n  Premature end to input file encountered : node index "  << " of file " << d_SpringInfoFile << std::endl);
        }
        else // read data into d_no_axial_node, d_no_circ_node, d_no_radial_node;
        {   d_vec_mesh_info.resize(3);
            line_string = discard_comments(line_string);
            std::istringstream line_stream(line_string);

            line_stream>>d_vec_mesh_info[0];
            line_stream>>d_vec_mesh_info[1];
            line_stream>>d_vec_mesh_info[2];
            pout << "++++++++++++++				 +++++++++++++++\n";
            pout <<"      ++++ Number of node in Axial  Direction = "<< d_vec_mesh_info[0] <<"\n "
                <<"      ++++ Number of node in Circ.  Direction = "<<d_vec_mesh_info[1]<<"\n "
               <<"      ++++ Number of node in Radial Direction = "<<d_vec_mesh_info[2]<<"\n";
        }

 	// >> added 02/28/2013 inputfile to the d_check_timeSteps
	if (!std::getline(file_stream, line_string))
        {
            d_check_timeSteps=0;// no dumping file
	   pout<< "+++++++++++++++++!! do not dump spring info! +++++++++++++++++\n";
        }
        else // read data into d_check_timeStep
        {   
            line_string = discard_comments(line_string);
            std::istringstream line_stream(line_string);
	    line_stream>>d_check_timeSteps;
            
            pout << "++++++++++++++				 +++++++++++++++\n";
            pout <<"      ++++ Dump Structure info. every steps= "<< d_check_timeSteps <<"\n ";
	    // added 03/09/2013 >> check the directory; we created a logfile: check.log 
	    
        }
	
	// << added 02/28/2013 add input>> d_check_timeStep

        // Close the input file.
	
        file_stream.close();

    }// TimeDependentSpring
    TimeDependentSpring::registerUpdateSpringFcn(0,&default_update_spring_fcn);
    return;
}

void
TimeDependentSpring::updateSpringConstants(
       // tbox::Pointer<hier::PatchHierarchy<NDIM> > const  hierarchy,
       // LDataManager* const lag_manager,
        const double current_time,
        const double dt,
	const int iteration_num)
{
    // >> added 02/27/2013 to deal with dumping springs into files
	bool do_check=false;
    if ((d_check_timeSteps>0) && (iteration_num%d_check_timeSteps== 0))
	 do_check=true;
    std::vector<DumpSpring> dump_springs;
    DumpSpring cur_dump_spring;
    // << added 02/27/2013
    const double tnext = current_time; //+dt-> in this case we will miss the t=0, first spring activation;
    
    int no_contracted_spring=0;
    int no_relaxed_spring=0;
    int no_shortening_spring=0;
    int no_lengthening_spring =0;  //
    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx =d_lag_manager->getLNodePatchDescriptorIndex();

    // find the spring that we need to modify
    const int finest_ln =d_hierarchy->getFinestLevelNumber();
    tbox::Pointer<hier::PatchLevel<NDIM> > level =d_hierarchy->getPatchLevel(finest_ln);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const  Pointer<LNodeSetData> lag_node_index_data  = patch->getPatchData(lag_node_index_idx);
        for(LNodeSetData::DataIterator it = lag_node_index_data->data_begin(patch_box);
            it != lag_node_index_data->data_end(); ++it)
        {
            LNode* const node_idx = *it;
            const int lag_idx = node_idx->getLagrangianIndex();

            // first check the node: whether the lag_idx belongs the lag node we need to modify, notice here this is the global index
            int spring_fcn_idx=-1; // this is the function index, -1 means not updated
            int spring_info_idx=IndexInfoNodeIn(lag_idx, spring_fcn_idx); // index of spring info that node is in
            if (spring_info_idx>-1) // -1, not in range;0 1,2 means the index of spring info
            {
                //pout<<"spring index is: "<< spring_fcn_idx<< "\n";
      IBAMR::IBSpringForceSpec* force_spec = node_idx->getNodeDataItem<IBAMR::IBSpringForceSpec>();              

                if (force_spec!= NULL) // we get the master node associate the active part,
                {   const int num_springs=force_spec->getNumberOfSprings();
                    const std::vector< int >  vec_fun_index = force_spec->getForceFunctionIndices ();
                    const int   mast_idx=force_spec->getMasterNodeIndex();
                    std::vector< double > & vec_rest_length=force_spec->getRestingLengths ();
                    std::vector< double> & vec_stiffnesses=force_spec->getStiffnesses ();
                    const std::vector< int > vec_slave_Nodes=force_spec->getSlaveNodeIndices ();
           /*         // check consistency of size of force_spec

            */
                    for (int it_springs=0;it_springs<num_springs;++it_springs)
                    {// second check: we check the force fcns, a unique marker for the active part
                        if ( vec_fun_index[it_springs]==spring_fcn_idx) //
                        {
                            // >> added 03/23/2013, get current state of spring, based on rest_length and bounds
                            int cur_state=GetCurrentState(spring_info_idx, vec_rest_length[it_springs]);
                            // << added 03/23/2013
                            int iflag=(d_update_spring_fcn_map[d_wave_info[spring_info_idx].idx_wave_algorithm])
                                    (d_wave_info[spring_info_idx].vec_wave_parameter,
                                     d_vec_mesh_info,
                                     mast_idx,
                                     vec_slave_Nodes[it_springs],
                                     vec_stiffnesses[it_springs],
                                     vec_rest_length[it_springs],tnext,cur_state);
                            //	no_actual_updated = no_actual_updated + iflag;
                            // if (Isupdate>0) no_actual_updated=++;
                            if (iflag==1) // contracted spring
                                no_contracted_spring++;
                            else if(iflag==2) // relaxed spring
                                no_relaxed_spring++;
                            else if(iflag==3)
                                no_shortening_spring++;
                            else if(iflag==4)
                                no_lengthening_spring++;
                            // >> added 02/27/2013 to store the spring info on this current processor;

                            if (do_check) // record and store spring info
                            {cur_dump_spring.id_struc=spring_info_idx;
                                cur_dump_spring.master_id=mast_idx;
                                cur_dump_spring.slave_id=vec_slave_Nodes[it_springs];
                                cur_dump_spring.rest_length=vec_rest_length[it_springs];
                                cur_dump_spring.spring_constant=vec_stiffnesses[it_springs];
                                cur_dump_spring.flag_id=iflag;
                                dump_springs.push_back(cur_dump_spring);
                            }
	        // << added 02/27/2013 to store the spring info.



				

                        }
                    }
                }
            }
        }
    }
   int sum_red1=SAMRAI::tbox::SAMRAI_MPI::sumReduction (no_contracted_spring);
   int sum_red2=SAMRAI::tbox::SAMRAI_MPI::sumReduction (no_relaxed_spring);
   int sum_red3=SAMRAI::tbox::SAMRAI_MPI::sumReduction (no_shortening_spring);
   int sum_red4=SAMRAI::tbox::SAMRAI_MPI::sumReduction (no_lengthening_spring);

   pout<<"[update_spring]: Number of contracted springs = " << sum_red1 << "\n";
   pout<<"[update_spring]: Number of relaxed springs= " << sum_red2 << "\n";
   pout<<"[update_spring]: Number of shortening springs= " << sum_red3 << "\n";
   pout<<"[update_spring]: Number of lengthening springs= " << sum_red4 << "\n";
   // >> added 
   if(do_check && dump_springs.size()>0)
{    // dir= ./CheckSpring 
	
    std::string file_name = "./CheckSpring/";
    char temp_buf[128];
    sprintf(temp_buf, "%06d.spring.proc%3d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    // begin to write to files
    std::ofstream DumpStream (file_name.c_str(),std::fstream::out);
    DumpStream << dump_springs.size() << " //struc_id, mast_id, slave_id, sp_const, rest_leng, iflag" << std::endl; 
    for (std::vector<DumpSpring>::size_type it=0; it<dump_springs.size(); it++)
    {DumpStream << dump_springs[it].id_struc << '\t' 
		<< dump_springs[it].master_id<< '\t'
		<< dump_springs[it].slave_id<< '\t'
		<< dump_springs[it].spring_constant<< '\t'
		<< dump_springs[it].rest_length<< '\t'
		<< dump_springs[it].flag_id << '\t'
		<< std::endl;
    }
    DumpStream.close();
    dump_springs.clear();
    pout<<"[update_spring]: dump to file..........\n";	
}
  
  // pout<< "[update spring]: done" << "\n";   
   return;  
  
}

int
TimeDependentSpring::IndexInfoNodeIn(const int node_idx, int& spring_fn)
{
    for (int k_range=0;k_range<d_no_spring_info;++k_range)
    {
        if ((node_idx>=d_spring_info[k_range].idx_min_node) && (node_idx<=d_spring_info[k_range].idx_max_node))
        {spring_fn= d_spring_info[k_range].idx_spring_fn;
            return k_range;
        }
    }


    return -1; // do not find the node

}// to check whether the lag. node is in the range of node with active spring

// to solve the bug: the change of rest_length may alter when regrid happens
// 1: contracted state, 0; normal state, -1: more relaxed state;
int 
TimeDependentSpring::GetCurrentState(int spring_info_idx, double rest_leng)
{
  if (rest_leng>d_bound_check[spring_info_idx].upper_check)
    return -1; // more relaxed, abnormal case
  else if(rest_leng<d_bound_check[spring_info_idx].lower_check)
    return 1; // contracted 
    
  return 0;
      
    
}
// 

}//TimeDependentSpring.C
