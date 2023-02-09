
#include "debug.hpp"

// bool DEBUG = false;
string dir = "Output/";
ofstream fileSummary(dir+"summary.txt");
ofstream fileModels(dir+"models.txt");
ofstream fileSolution(dir+"data.txt");

string indent = "";

string now() {
    /*
     *
     */

    time_t rawTime;
    struct tm* timeInfo;
    time(&rawTime);
    timeInfo = localtime(&rawTime);
    char timeStamp [80];
    strftime(timeStamp, sizeof (timeStamp), "%H:%M:%S", timeInfo);
    return string(timeStamp);
}

void dbg(ostream &out, string s, bool newline, bool time) {
    /*
     *
     */

    //if(!DEBUG) return;

    if(time) out << now() << " ";
    out << indent << s;
    if(newline) out << "\n";
}

void dbgc(ostream &out, string s) {
    /*
     *
     */

    //if(!DEBUG) return;
    out << s;
}

void dbgnl(ostream &out, string s) {
    /*
     *
     */

    //if(!DEBUG) return;
    dbgc(out,s+"\n");
}

void in() {
    /*
     *
     */

    //if(!DEBUG) return;
    indent = indent + "  ";
}

void out() {
    /*
     *
     */

    //if(!DEBUG) return;
    indent = indent.substr(0, indent.size()-2);
}

string ind() {
    /*
     *
     */

    //if(!DEBUG) return "";
    return indent;
}