//create a directory
#ifndef CREATEDDIRECTORY_H_
#define CREATEDDIRECTORY_H_
#include "clibrary.h"
#include <ctype.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

namespace createdirectory
{
    int mkpath(string s,mode_t mode=0755)
    {
        size_t pre=0,pos;
        std::string dir;
        int mdret,errno_0;
        
        if(s[s.size()-1]!='/')
        {
            // force trailing / so we can handle everything in loop
            s+='/';
        }
        
        errno_0 = errno;
        while((pos=s.find_first_of('/',pre))!=string::npos)
        {
            errno = errno_0;
            
            dir=s.substr(0,pos++);
            pre=pos;
            if(dir.size()==0) continue; // if leading / first time is 0 length
            #if defined(WIN32)
            if((mdret=::mkdir(dir.c_str())) && errno!=EEXIST) // && pos!=(s.size()-2)
            {
                return mdret;
            }
            #else
            if((mdret=::mkdir(dir.c_str(),mode)) && errno!=EEXIST) // && pos!=(s.size()-2)
            {
                return mdret;
            }
            #endif
            if((pos==s.size()) && errno==EEXIST) mdret = 0;
            
        }
        return mdret;
    }
}
#endif //CREATEDDIRECTORY_H_
