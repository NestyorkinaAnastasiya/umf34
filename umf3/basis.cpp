/*basis.cpp*/
#include "basis.h"

namespace basis
{
	Basis::Basis()
	{
		array <function<double(double)>, nFunc1D> phi_;
		array <function<double(double)>, nFunc1D> dphi_ksi;
		phi_[0] = [](double ksi) { return 1 - ksi; };
		phi_[1] = [](double ksi) { return ksi; };
		dphi_ksi[0] = [](double ksi) { return -1; };
		dphi_ksi[1] = [](double ksi) { return  1; };
		//указатели на функции вычисления базисных функций в точке
		phi[0] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[0](etta); };
		phi[1] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[0](etta); };
		phi[2] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[1](etta); };
		phi[3] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[1](etta); };
		//указатели на функции вычисления d/dksi базисных функций в точке
		dphiksi[0] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[0](etta); };
		dphiksi[1] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[0](etta); };
		dphiksi[2] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[1](etta); };
		dphiksi[3] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[1](etta); };
		//указатели на функции вычисления d/detta базисных функций в точке
		dphietta[0] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[0](etta); };
		dphietta[1] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[0](etta); };
		dphietta[2] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[1](etta); };
		dphietta[3] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[1](etta); };
	}
}
		
		/*
		array <function<double(double)>, nFunc1D> phi_;
		phi_[0] = [](double ksi) { return 2 * (ksi - 0.5) * (ksi - 1); };
		phi_[1] = [](double ksi) { return -4 * ksi * (ksi - 1); };
		phi_[2] = [](double ksi) { return 2 * ksi * (ksi - 0.5); };

		array <function<double(double)>, nFunc1D> dphi_ksi;
		dphi_ksi[0] = [](double ksi) { return 4 * ksi - 3; };
		dphi_ksi[1] = [](double ksi) { return  -8 * ksi + 4; };
		dphi_ksi[2] = [](double ksi) { return 4 * ksi - 1; };

		phi[0] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[0](etta); };
		phi[1] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[0](etta); };
		phi[2] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[0](etta); };
		phi[3] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[1](etta); };
		phi[4] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[1](etta); };
		phi[5] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[1](etta); };
		phi[6] = [phi_](double ksi, double etta) { return phi_[0](ksi) * phi_[2](etta); };
		phi[7] = [phi_](double ksi, double etta) { return phi_[1](ksi) * phi_[2](etta); };
		phi[8] = [phi_](double ksi, double etta) { return phi_[2](ksi) * phi_[2](etta); };

		dphiksi[0] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[0](etta); };
		dphiksi[1] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[0](etta); };
		dphiksi[2] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[0](etta); };
		dphiksi[3] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[1](etta); };
		dphiksi[4] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[1](etta); };
		dphiksi[5] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[1](etta); };
		dphiksi[6] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[0](ksi) * phi_[2](etta); };
		dphiksi[7] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[1](ksi) * phi_[2](etta); };
		dphiksi[8] = [phi_, dphi_ksi](double ksi, double etta) { return dphi_ksi[2](ksi) * phi_[2](etta); };

		dphietta[0] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[0](etta); };
		dphietta[1] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[0](etta); };
		dphietta[2] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[0](etta); };
		dphietta[3] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[1](etta); };
		dphietta[4] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[1](etta); };
		dphietta[5] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[1](etta); };
		dphietta[6] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[0](ksi) * dphi_ksi[2](etta); };
		dphietta[7] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[1](ksi) * dphi_ksi[2](etta); };
		dphietta[8] = [phi_, dphi_ksi](double ksi, double etta) { return phi_[2](ksi) * dphi_ksi[2](etta); };*/
	
