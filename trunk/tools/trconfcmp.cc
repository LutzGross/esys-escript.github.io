

#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <set>
#include <algorithm>
#include <vector>
#include <cctype>

using namespace std;

class rec{
public:    
    string name;
    string type;
    string value;
    bool operator==(const rec& r1) const {
        return r1.name==this->name;
    }
    bool operator<(const rec& r1) const {
        return r1.name<this->name;
    }
    
    friend ostream& operator<<(ostream& os, const rec& r);

};

ostream& operator<<(ostream& os, const rec& r) {
    os << '[' << r.name << '|' << r.type << "](" << r.value << ")" << endl;
    return os;    
}



auto loadmap(ifstream& ifs) {
    string line;
    set<rec> res;
    while(getline(ifs, line), ifs.good()) {
        auto pos=line.find("-D");
        if (pos==string::npos) {
            continue;
        }
        while (pos+2<line.size() && (line[pos+2]==' ')) {
            pos++;
        }
        pos+=2;
        auto tpos=line.find(':', pos);
        auto eqpos=line.find('=', pos);
        auto end=line.find_last_of('\\');
        if (eqpos==string::npos) {
            continue;
        }
        if (end==string::npos) {
            end=line.size()-1;
        } else {
            end-=1; // move back from the backslash
        }
        while (isspace(line[end])) {
            end--;
        }
        rec r;
        r.value=line.substr(eqpos+1, end-eqpos);
        if (tpos!=string::npos) {
            r.type=line.substr(tpos+1, eqpos-tpos-1);
        } else {
            tpos=eqpos;     // type is missing so we pretend = and : coincide
        }
        r.name=line.substr(pos,tpos-pos);
        res.insert(r);
        //cout << r.name << ' ' << r.type << " " << r.value << endl;
    }
    return res;
}


int main(int argc, char** argv) {
    if (argc!=3) {
        cerr << argv[0] << " file1 file2\n";
        return 1;
    }
    string line;
    ifstream f1(argv[1]);
    ifstream f2(argv[2]);
    if (!f1 || !f2) {
        cerr << "Can't open file\n";
        return 2;
    }
    set<rec> r1=loadmap(f1);
    set<rec> r2=loadmap(f2);

    vector<rec> all;
    vector<rec> only_left;
    vector<rec> only_right;
    vector<pair<rec,rec>> diff_types;
    vector<pair<rec,rec>> diff_values;
    all.reserve(r1.size()+r2.size());
    set_union(r1.begin(), r1.end(), r2.begin(), r2.end(), inserter(all, all.begin()));
    
    for (auto it=all.begin();it!=all.end();++it) {
        const rec* lhs=0;
        const rec* rhs=0;
        for (auto il=r1.begin();il!=r1.end();++il) {
            if (il->name==it->name) {
                lhs=&(*il);
            }
        }
        for (auto il=r2.begin();il!=r2.end();++il) {
            if (il->name==it->name) {
                rhs=&(*il);
            }
        }
        if (rhs==0) {
            only_left.push_back(*lhs);
        } else if (lhs==0) {
            only_right.push_back(*rhs);
        } else if (lhs->value!=rhs->value) {
            diff_values.push_back(pair<rec,rec>(*lhs, *rhs));
        } else if (lhs->type!=rhs->type) {
            diff_types.push_back(pair<rec,rec>(*lhs, *rhs));            
        } else {
        }
    }
    if (!only_left.empty())
    {
        cout << "Only in LHS:\n";
        for (auto it=only_left.begin();it!=only_left.end();++it) {
            cout << "    -D" << it->name;
            if (!it->type.empty()) {
                cout << ':' << it->type;
            }
            cout << '=' << it->value << " \\\n";
        }
    }
    if (!only_right.empty())
    {
        cout << "Only in RHS:\n";
        for (auto it=only_right.begin();it!=only_right.end();++it) {
            cout << "    -D" << it->name;
            if (!it->type.empty()) {
                cout << ':' << '{' << it->type << '}';
            }
            cout << '=' << it->value << " \\\n";
        }
    }
    if (!diff_values.empty())
    {
        cout << "Different values:\n";
        for (auto it=diff_values.begin();it!=diff_values.end();++it) {
            cout << "    " << it->first.name;
            if (!it->first.type.empty()) {
                cout << ':' << it->first.type;
            } else if (!it->second.type.empty()) {
                cout << ':' << it->second.type;                
            }
            cout << " LHS=" << it->first.value << " RHS=" << it->second.value << endl;
        }
    }
    if (!diff_types.empty()) {
        cout << "Different types:\n";
        for (auto it=diff_types.begin();it!=diff_types.end();++it) {
            cout << "    " << it->first.name;
            cout << " LHS=" << it->first.type << " RHS=" << it->second.type << endl;
        }
    }
}
