#include <iostream>
#include <string>
#include <array>
#include <regex>
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"

class MediaPlanInfo 
{
    public:
        // construct:
        MediaPlanInfo(const std::string mp_info_str);
        
        // members:
        std::string orig_str;
        unsigned int id;
        unsigned int inventory_id;
        std::string day_part;
        unsigned int day_part_num;
        unsigned int year;
        unsigned int month;
        unsigned int day;
        unsigned int hour;
        unsigned int min;
        unsigned int sec;
        unsigned int creative_id;
        unsigned int creative_duration;
        unsigned int var_type; 
};

// non-members:
std::ostream& operator << (std::ostream& out, const MediaPlanInfo& mp_info) 
{
    out << "MediaPlanInfo{"
        << "\n    orig_str          : '" << mp_info.orig_str << "'"
        << "\n    id                : "  << mp_info.id
        << "\n    inventory_id      : "  << mp_info.inventory_id
        << "\n    day_part          : "  << mp_info.day_part
        << "\n    day_part_num      : "  << mp_info.day_part_num
        << "\n    year              : "  << mp_info.year
        << "\n    month             : "  << mp_info.month
        << "\n    day               : "  << mp_info.day
        << "\n    hour              : "  << mp_info.hour
        << "\n    min               : "  << mp_info.min
        << "\n    sec               : "  << mp_info.sec
        << "\n    creative_id       : "  << mp_info.creative_id
        << "\n    creative_duration : "  << mp_info.creative_duration
        << "\n    var_type          : "  << mp_info.var_type
        << "\n}";
    return out;
}

// class implementation (to put in cpp file eventually):
MediaPlanInfo::MediaPlanInfo(const std::string mp_info_str)
    : orig_str(mp_info_str) 
    , id(9999)
    , inventory_id(9999)
    , day_part("bogus")
    , day_part_num(9999)
    , year(2014)
    , month(9999)
    , day(9999)
    , hour(9999)
    , min(9999)
    , creative_id(9999)
    , creative_duration(9999)
    , var_type(9999)
{
    std::regex re(R"((\d{3})-(\d{4})\[(\w+)\((\d)\)\]@(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})\((\d{2})-(\d{1,3})s\)#(\d))");
    std::smatch match;
    std::regex_search(orig_str, match, re);

    // check boundary conditions:
    if (match.empty())
    {
        throw std::runtime_error("[MediaPlanInfo] Error : could not match media plan string"); 
    } 

    // fill attributes:
    try
    {
        id                = boost::lexical_cast<unsigned int>(match[1 ].str());
        inventory_id      = boost::lexical_cast<unsigned int>(match[2 ].str());
        day_part          = match[3 ].str() ;
        day_part_num      = boost::lexical_cast<unsigned int>(match[4 ].str());
        year              = 2014; // hard coded in field?
        month             = boost::lexical_cast<unsigned int>(match[5 ].str());
        day               = boost::lexical_cast<unsigned int>(match[6 ].str());
        hour              = boost::lexical_cast<unsigned int>(match[7 ].str());
        min               = boost::lexical_cast<unsigned int>(match[8 ].str());
        sec               = boost::lexical_cast<unsigned int>(match[9 ].str());
        creative_id       = boost::lexical_cast<unsigned int>(match[10].str());
        creative_duration = boost::lexical_cast<unsigned int>(match[11].str());
        var_type          = boost::lexical_cast<unsigned int>(match[12].str());
    }
    catch (boost::bad_lexical_cast& e)
    {
        // better error messege
        std::cerr << "[MediaPlanInfo] Error : lexical cast failed in constructor\n"; 
        throw e;
    }
}

int main()
try
{
    const std::array<std::string, 3> strs = {{
        "165-2964[Morning(3)]@0915080751(61-5s)#0", 
        "165-2964[Morning(3)]@0915080751(62-15s)#0", 
        "165-2961[Early_Morning(2)]@0915063000(61-5s)#1"
    }};
    for (const std::string& str : strs)
    {
        MediaPlanInfo mp_info(str);
        std::cout << "mp_info = " << mp_info << std::endl; 
    }
    return 0;
}
catch (std::exception& e)
{
    // Syntax error in the regular expression
    std::cout << "[main error]" << e.what() << std::endl;
}

