
/*****************************************************************************
*
* Copyright (c) 2017-2018 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/



#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <boost/lexical_cast.hpp>

/* From http://www.unidata.ucar.edu/software/netcdf/docs/file_format_specifications.html
The grammar for the header is:

header       = magic  numrecs  dim_list  gatt_list  var_list
= CDF? followed by 4 bytes(which we don't care about)
We can just search for the relevant tags:


*/

union esV
{
    unsigned char b;
    short s;
    int i;
    float f;
    double d;
};


class esncdfValues
{
public:
    char type;
    std::string svalue;
    std::vector<union esV> vec;  
};

namespace
{
    char NC_DIM_TAG=10;
    char NC_VAR_TAG=11;
    char NC_ATT_TAG=12;
    
    
    
    enum nc_type
    {
      NC_BYTE=1,
      NC_CHAR=2,
      NC_SHORT=3,
      NC_INT=4,
      NC_FLOAT=5,
      NC_DOUBLE=6
    };
    
    int offset;
    
char getChar(std::ifstream& ifs)
{
    offset++;
    return ifs.get(); 
}




typedef std::vector<std::pair<std::string, esncdfValues> >  vsu; 

}


using namespace std;
using namespace boost;








// Holds a list of name value pairs
class esncdfCollection
{
public:
    void add(const std::string& k, const esncdfValues& v)
    {
        vec.push_back(std::pair<std::string, esncdfValues>(k,v));
    }
    void add_int(const std::string& k, int i)
    {
	esncdfValues v;
	v.type=NC_INT;
	union esV es;
	es.i=i;
	v.vec.push_back(es);
        vec.push_back(std::pair<std::string, esncdfValues>(k,v));
    }
    
/*    void add(const std::string& k, const std::string& s)
    {
	esncdfValues u;
	u.type=NC_CHAR;
	u.svalue=s;
        vec.push_back(std::pair<std::string, esncdfValues>(k,u));
    }  */  
    vsu vec;
    // Here would be extra information regarding variables
    std::string name;       
    std::string type;
};


class esncdfVariables
{
public:    
    std::vector<esncdfCollection> vars;    
    void add(const esncdfCollection& c)
    {
        vars.push_back(c);
    }
};


class esncdfHeader
{
public:    
    esncdfCollection dimensions;
    esncdfCollection attributes;     // global attributes
    esncdfVariables variables;       // including specific attributes
};


bool read_list(std::ifstream& ifs, char t, bool offset32, esncdfCollection& coll);


bool gobble(ifstream& ifs, int s)
{
    for (int i=0;i<s;++i)
    {
        getChar(ifs);
    }
    return ifs.good();
}


bool read_int(std::ifstream& ifs, int& res)
{
    char buff32[4];
    for (int i=0;i<4;++i)
    {
        buff32[i]=getChar(ifs);
    }
    if (!ifs)
    {
        return false;
    }
    res=0;
    for (int i=0;i<4;++i)
    {
        res=(res << 8);
        res+=(unsigned char)buff32[i];
    }
    return true;
}

bool read_item_header(std::ifstream& ifs, std::string& iname)
{
    // format for name is 32 bit length followed by chars (padded to 4byte boundary)
    int len=0;
    if (!read_int(ifs,len))
    {
        return false;
    }
    for (int i=0;i<len;++i)
    {
        char c=getChar(ifs);
        if (c!=0)       // if anyone puts \x00 in their string
        {               // ...
            iname+=c;
        }
    }
    // now we need to skip padding
    if (len%4!=0)
    {
	for (int i=0;i<4-len%4;++i)
	{
	    char c=getChar(ifs);
	}
    }
    if (!ifs)
    {
        return false;
    }
    return true;    
}

bool read_item_dim(std::ifstream& ifs, std::string& iname, esncdfCollection& dims)
{
    if (!read_item_header(ifs, iname))
    {
        return false;
    }
    int len=0;
    if (!read_int(ifs, len))
    {
        return false;
    }
    dims.add_int(iname, len);
    return true;
}

bool read_tag(std::ifstream& ifs, char& t)
{
    for (int i=0;i<3;++i)
    {
        t=getChar(ifs); 
	if (t!=0)
	{
	    return false;
	}
    }
    t=getChar(ifs);
    return ifs.good();  
}

esncdfValues get_string_value(char tag, int len, char* bytes)
{
    esncdfValues res;
    res.type=tag;
    switch (tag)
    {
        case NC_CHAR:  res.svalue=std::string(bytes); break;
        case NC_BYTE:
	{
	    union esV e;
            for (int i=0;i<len;++i)
            {
		e.b=bytes[i];
		res.vec.push_back(e);
            }
            break;
	}
        case NC_SHORT:
        {
	    union esV e;
	    for (int i=0;i<len;++i)
	    {
		short v=(bytes[0] << 8)+bytes[1];
	        e.s=v;
		res.vec.push_back(e);
	    }
            break;
        }
        case NC_INT:
        {
	    union esV e;
	    for (int i=0;i<len;++i)
	    {
		int v=(bytes[0] << 24) +(bytes[1] << 16) + (bytes[2] << 8) +bytes[3];	      
	        e.i=v;
		res.vec.push_back(e);
	    }
            break;
        }
        case NC_FLOAT:      // here we assume the correct floating point std
        {
	    union esV e;
	    float f=-0.0;
	    if (reinterpret_cast<unsigned char*>(&f)[0]!=0)	// big endian (yes we assume IEEE floating point
	    {
		for (int i=0;i<len;++i)
		{
		    float f=0;		// assumes floats are 32bit!
		    unsigned char* v=reinterpret_cast<unsigned char*>(&f);
		    for (int j=0;j<4;++j)
		    {
			v[j]=bytes[4*i+j];
		    }
		    e.f=f;
		    res.vec.push_back(e);
		}
	    }
	    else
	    {
		for (int i=0;i<len;++i)
		{
		    float f=0;		// assumes floats are 32bit!
		    unsigned char* v=reinterpret_cast<unsigned char*>(&f);
		    for (int j=0;j<4;++j)
		    {
			v[j]=bytes[4*i+3-j];
		    }
		    e.f=f;
		    res.vec.push_back(e);
		}
	    }
	    break;
        }
        case NC_DOUBLE:
        {
	    union esV e;
	    double f=-0.0;
	    if (reinterpret_cast<unsigned char*>(&f)[0]!=0)	// big endian (yes we assume IEEE floating point
	    {
		for (int i=0;i<len;++i)
		{
		    double f=0;
		    unsigned char* v=reinterpret_cast<unsigned char*>(&f);
		    for (int j=0;j<8;++j)
		    {
			v[j]=bytes[8*i+j];
		    }
		    e.d=f;
		    res.vec.push_back(e);
		}
	    }
	    else
	    {
		for (int i=0;i<len;++i)
		{
		    double f=0;
		    unsigned char* v=reinterpret_cast<unsigned char*>(&f);
		    for (int j=0;j<8;++j)
		    {
			v[j]=bytes[8*i+7-j];
		    }
		    e.d=f;
		    res.vec.push_back(e);
		}
	    }	  
	    break;	  
        }
        
        
        default:
	    res.type=NC_CHAR;
	    res.svalue="?????";
    }    
    return res;  
}

// structure for an attribute is:
//  name  nc_type  nelems  [values ...]
bool read_item_attr(std::ifstream& ifs, std::string& iname, esncdfCollection& coll)
{
    if (!read_item_header(ifs, iname))
    {
        return false;
    }    
    int len=0;
    char tag;
    if (!read_tag(ifs, tag))
    {
        return false;
    }
    if ((tag<NC_BYTE) || (tag>NC_DOUBLE))	// invalid type
    {
        return false;
    }
    if (!read_int(ifs, len))
    {
        return false;
    }
    // for now I'm just going to skip over the values
    // to do that, I need to know how big each type is
    int size=0;
    switch (tag)
    {
        case NC_SHORT: size=2; break;
        case NC_INT: size=4; break;
        case NC_FLOAT: size=4; break;
        case NC_DOUBLE: size=8; break;
        default:
        size=1;
    }
    
    int skip=len*size;
    // now we need to skip padding
    if (len*size%4!=0)
    {
        skip+=4-len*size%4;
    }
    char buff[skip+1];
    buff[skip]=0;
    for (int i=0;i<skip;++i)
    {
        buff[i]=getChar(ifs);
    }
    if (!ifs.good())
    {
        return false;
    }
    esncdfValues vs=get_string_value(tag, len, buff);
    coll.add(iname, vs);    
    return true;
}





// structure for a variable is:
// name  nelems  [dimid ...]  vatt_list  nc_type  vsize  begin
bool read_item_var(std::ifstream& ifs, std::string& iname, bool offset32, esncdfVariables& vars)
{
    if (!read_item_header(ifs, iname))
    {
        return false;
    }    
    int dim=0;
    if (!read_int(ifs, dim))
    {
        return false;
    }
    for (int i=0;i<dim;++i)	// just disgarding these at the moment
    {
      int v;
      if (!read_int(ifs, v))
      {
        return false;
      } 
    }
    esncdfCollection current;
    current.name=iname;
    if (!read_list(ifs, NC_ATT_TAG, offset32, current))
    {
        return false;
    }
    // now determine the type
    char tag;
    if (!read_tag(ifs, tag))
    {
        return false;
    }
    if ((tag<NC_BYTE) || (tag>NC_DOUBLE))	// invalid type
    {
        return false;
    }
    switch(tag)
    {
    case NC_CHAR: current.type="char"; break;
    case NC_BYTE: current.type="byte";break;
    case NC_SHORT: current.type="short";break;
    case NC_INT: current.type="int";break;
    case NC_FLOAT: current.type="float";break;
    case NC_DOUBLE: current.type="double";break;
    default:
        current.type="????";
    }
    int len;
    if (!read_int(ifs, len))
    {
        return false;
    }
    // now we need to get the offset of the variable in the file
    // which we will promptly discard
    int gobcount=offset32?4:8;
    if (!gobble(ifs, gobcount))
    {
        return false;
    }
    vars.add(current);
    return true;
}

// lists<T> are defined as: ABSENT | NC_DIMENSION  nelems  [T ...]
// where T in {dim, attr, var}
bool read_list_vars(std::ifstream& ifs, bool offset32, esncdfVariables& coll)
{
    char b;
    // we expect either 8 zeros indicating no list or the expected tag followed
    // by list information
    for (int i=0;i<3;++i)
    {
        
        if (!ifs || ((b=getChar(ifs))!=0))
        {
            return false;
        }
    }
    b=getChar(ifs);
    if (!ifs || (b!=0 && b!=NC_VAR_TAG)) // didn't get tag we expected
    {
        return false;
    }
    if (b==0)   // expecting 4 more zeros
    {
        for (int i=0;i<4;++i)
        {
            b=getChar(ifs);
            if (!ifs || b!=0)
            {
                return false;
            }
        }
        return true;		// This list is absent
    }
    // we must have seen the correct tag
    // first thing will be number of items
    int count=0;
    if (!read_int(ifs, count))
    {
        return false;
    }
    for (int i=0;i<count;++i)
    {
        std::string itemname;
        if (!read_item_var(ifs, itemname, offset32, coll))
        {
            return false;
        }
    }
    return true;
}




// lists<T> are defined as: ABSENT | NC_DIMENSION  nelems  [T ...]
// where T in {dim, attr, var}
bool read_list(std::ifstream& ifs, char t, bool offset32, esncdfCollection& coll)
{
    char b;
    // we expect either 8 zeros indicating no list or the expected tag followed
    // by list information
    for (int i=0;i<3;++i)
    {
        
        if (!ifs || ((b=getChar(ifs))!=0))
        {
            return false;
        }
    }
    b=getChar(ifs);
    if (!ifs || (b!=0 && b!=t)) // didn't get tag we expected
    {
        return false;
    }
    if (b==0)   // expecting 4 more zeros
    {
        for (int i=0;i<4;++i)
        {
            b=getChar(ifs);
            if (!ifs || b!=0)
            {
                return false;
            }
        }
        return true;		// This list is absent
    }
    // we must have seen the correct tag
    // first thing will be number of items
    int count=0;
    if (!read_int(ifs, count))
    {
        return false;
    }
    if (t==NC_DIM_TAG) 
    {
        for (int i=0;i<count;++i)
        {
            std::string itemname;
            if (!read_item_dim(ifs, itemname, coll))
            {
                return false;
            }
        }
    }
    else if (t==NC_ATT_TAG)
    {
        for (int i=0;i<count;++i)
        {
            std::string itemname;
            if (!read_item_attr(ifs, itemname, coll))
            {
                return false;
            }
        }
    }
    else
    {
        return false;
    }
    return true;
}

bool read_header(const char* fname, esncdfHeader& header)
{
    ifstream ifs(fname);
    if (!ifs)
    {
        return false;
    }
    gobble(ifs,3);	// skip over "CDF"
    bool offset32=(getChar(ifs)==1);
    gobble(ifs, 4);
    if (!ifs)
    {
        return false;
    }
    if (!read_list(ifs, NC_DIM_TAG, offset32, header.dimensions))        
    {
        return false;
    }
    cout << "Dimensions:\n";
    cout << "Attributes:\n";    
    if (!read_list(ifs, NC_ATT_TAG, offset32, header.attributes))        
    {
        return false;
    }

    cout << "Variables:\n";    
    if (!read_list_vars(ifs, offset32, header.variables)) 
    {
        return false;
    }   
    return true;
}

void dump_Values(std::ostream& os, const esncdfValues& u)
{
    if (u.type==NC_CHAR)
    {
        os << u.svalue;
    }
    else
    {
        for (int i=0;i<u.vec.size();++i)
	{
	    switch (u.type)
	    {
	    case NC_BYTE: 
	    {
		unsigned char x=u.vec[i].b;
	        os << hex << (unsigned short)(x) << dec << " "; break;
	    }
	    case NC_SHORT: os << u.vec[i].s << " "; break;
	    case NC_INT: os << u.vec[i].i << " "; break;
	    case NC_FLOAT: os << u.vec[i].f << " "; break;
	    case NC_DOUBLE: os << u.vec[i].d << " "; break;
	    default:
	      os << "????";
	    }
	}
    }
  
}

void print_header(const esncdfHeader& h)
{
    cout << "Dimensions:\n";
    for (auto it=h.dimensions.vec.begin();it!=h.dimensions.vec.end();++it)
    {
       const esncdfValues& v=it->second;
       cout << it->first << " = ";
       dump_Values(cout, v);
       cout << endl;
    }
    cout << "Attributes:\n";
    for (auto it=h.attributes.vec.begin();it!=h.attributes.vec.end();++it)
    {
       const esncdfValues& v=it->second;
       cout << it->first << " = ";
       dump_Values(cout, v);
       cout << endl;
    }
    cout << "Variables:\n";
    for (auto it=h.variables.vars.begin();it!=h.variables.vars.end();++it)
    {
       cout << it->type << " " << it->name << ":\n";
       for (auto jt=it->vec.begin();jt!=it->vec.end();++jt)
       {
	  cout << "  " << jt->first << " = ";
	  const esncdfValues& v=jt->second;
	  dump_Values(cout, v);
	  cout << endl;
       }
    }
}

int main(int argc, char** argv)
{
    offset=-1;
    if (argc==1)
    {
        std::cerr << "Usage: " << argv[0] << " file.nc\n";
        return 1;
    }
    esncdfHeader header;
    if (!read_header(argv[1], header))
    {
	std::cerr << "Failed\n";
    }
    cout << "-------------\n";
    print_header(header);
    return 0;

}
