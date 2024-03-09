#include <bits/stdc++.h>
using namespace std;
int main() {
  string S,s;
  cin>>S;
reverse(S.begin(),S.end());
int a=0;
bool b=true;
while (a<S.size()&&b)
{
  s="";
  for (int i = 0; i < 5; i++)
  {
    s=s+S.at(a+i);
  }
  if(s=="maerd"||s=="esare"){
    a+=5;
  }
  else{
    s=s+S.at(a+5);
    if(s=="resare"){
      a+=6;
    }
    else{
      s=s+S.at(a+6);
      if(s=="remaerd"){
      a+=7;
      }
      else{
        cout<<"NO"<<endl;
        b=false;
        break;
      }
    }
  }
}
if(b){
  cout<<"YES"<<endl;
}
return 0;
}