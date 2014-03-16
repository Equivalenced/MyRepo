#include "Squid2.h"

using namespace std;


int main(void) {
	Squid squid;
	squid.init0();
	squid.readIn0();
	squid.DoFFT();
	//squid.printFFT();
	squid.psd0();
	squid.doMedian();
	cout << "Squid Out - ";
	squid.write0();
	cout << "Complete" << endl;
	//squid.readIn();
	//squid.printData();
	//squid.write();
	//cout << "Press Enter to Continue ..." << endl;
	//cin.get();
	return (0);
}
