
int windows_size = 60;
int buffer;
bool odd;

void start_print(TString s){
    if( s.Length() % 2 ) odd = true;
    else odd = false;

    if(odd) s = " "+s + " ... Start! ";
    else s = " "+s + " ... Start!! ";

    buffer = ( windows_size-s.Length() )/2;
    cout<<endl;
	cout<<std::string(windows_size, '*')<<endl;
	cout<<std::string(buffer, '*')<<s<<std::string(buffer, '*')<<endl;
}

void line_print(){
    cout<<std::string(windows_size, '*')<<endl;
}


void running_print(TString s){
    if( s.Length() % 2 ) odd = true;
    else odd = false;

    if(odd) s = " "+ s;
    else s = " "+s + " ";

    buffer = ( windows_size-s.Length() )/2;
	cout<<"*"<<s<<std::string(2*buffer-2, ' ')<<"*"<<endl;
}

void finish_print(TString s){
    if( s.Length() % 2 ) odd = true;
    else odd = false;

    if(odd) s = " "+s + " ... Done!! ";
    else s = " "+s + " ... Done!!! ";

    buffer = ( windows_size-s.Length() )/2;
	cout<<std::string(buffer, '*')<<s<<std::string(buffer, '*')<<endl;
	cout<<std::string(windows_size, '*')<<endl;
}