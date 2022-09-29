#ifndef __TTM_FUNC__
#define __TTM_FUNC__ "pulse"

// all variables needed for this 2tm solution
namespace ttm_n
{
	double delay = 0.0;
	double tau = 0.0;
	double power = 0.0;
	
	double delta_sample = 100; // in number of layers
	double gamma_sample = 135.0;
	double Debye_sample = 30.0 * SI_milliElectronvolt;
	double c_ph_sample = 3.0e6;
    double g_ep_sample = 1.0e17;
	double2tm kappa_sample = {10.0, 1.0};
	
	double delta_substrate = 100; // in number of layers
	double gamma_substrate = 250.0;
	double Debye_substrate = 40.0 * SI_milliElectronvolt;
    double c_ph_substrate = 4.0e6;
    double g_ep_substrate = 1.0e17;
	double2tm kappa_substrate = {50.0, 1.0};
	
	double2tm kappa_interface = {1.0e9, 1.0e8};

    double2tm *predictor = nullptr; 
    double2tm *corrector = nullptr; 
    double2tm *result = nullptr;

	double funcEinstein(double kT, double Debye);
    double funcLaserPower(double time);
    int getStepSize(ttm_t ttm, double time, double deltaTime);
};

// make one time step with the Heun method
int ttm_t::update(double deltaTime)
{
    int status = STATUS_SUCCESS;

    const int width = this->size;

    // tmp variables
    double2tm qInterface = { 0 };
    double2tm qBulk = { 0 };
    double c_el, c_ph, g_ep, gradKappa, kappa, attenuation;
	
    // either predictor or corrector
	double2tm *activeField;

    // alias for better readability
    double2tm *kT_loc = this->kT_field_64bit;

    // sometimes smaller time steps are required for 2TM
    int stepMultiple = ttm_n::getStepSize(*this, this->time, deltaTime);
	deltaTime /= (double)stepMultiple;
	
	for(int step = 0; step < stepMultiple; step++)
    {
		memcpy(ttm_n::result, kT_loc, this->size * sizeof(double2tm));
		activeField = ttm_n::predictor;
		
        // predictor - corrector loop (0:predictor; 1:corrector)
        for (int integratorStep = 0; integratorStep < 2; integratorStep++)
        {
			memset(activeField, 0, this->size * sizeof(double2tm));
			
            // no flow condition
            kT_loc[0] = kT_loc[2];
			kT_loc[width - 1] = kT_loc[width - 3];

            // sample loop
			for(int z = 1; z < this->dim + 1; z++)
            {
				
				////////////////
				// ELECTRONS
				////////////////
				
                c_el = ttm_n::gamma_sample * kT_loc[z].e / SI_kBoltzmann;
				c_ph = ttm_n::c_ph_sample * ttm_n::funcEinstein(kT_loc[z].p, ttm_n::Debye_sample);
				g_ep = ttm_n::g_ep_sample;
                attenuation = exp(-(double)z / ttm_n::delta_sample);

				activeField[z].e += -g_ep * (kT_loc[z].e - kT_loc[z].p) * (deltaTime / c_el);
				activeField[z].e += ttm_n::funcLaserPower(time) * attenuation * SI_kBoltzmann * (deltaTime / c_el);
				activeField[z].e += -this->deltaQ_perStep.e * (SI_kBoltzmann / c_el) / stepMultiple;

                kappa = ttm_n::kappa_sample.e * (kT_loc[z].e / kT_loc[z].p);

                // special treatment for non-local terms at interface
				if (z == this->dim)
                {
                    qInterface.e = -ttm_n::kappa_interface.e  * (kT_loc[z + 1].e - kT_loc[z].e);
                    qBulk.e = -kappa * (kT_loc[z].e - kT_loc[z - 1].e);
					activeField[z].e += -(qInterface.e - qBulk.e) * (deltaTime / c_el);
				}
				else
                {
					activeField[z].e += kappa * (kT_loc[z + 1].e - 2.0 * kT_loc[z].e + kT_loc[z - 1].e) * (deltaTime / c_el);
					gradKappa = 0.5 * ttm_n::kappa_sample.e * (kT_loc[z + 1].e / kT_loc[z + 1].p - kT_loc[z - 1].e / kT_loc[z - 1].p);
					activeField[z].e += 0.5 * gradKappa * (kT_loc[z + 1].e - kT_loc[z - 1].e) * (deltaTime / c_el);
				}
				
				//////////////
				// PHONONS
				//////////////

				activeField[z].p += g_ep * (kT_loc[z].e - kT_loc[z].p) * (deltaTime / c_ph);
                activeField[z].p += -this->deltaQ_perStep.p * (SI_kBoltzmann / c_ph) / stepMultiple;

                // special treatment for non-local terms at interface
				if (z == this->dim)
                {
					qInterface.p = -ttm_n::kappa_interface.p * (kT_loc[z + 1].p - kT_loc[z].p);
					qBulk.p = -ttm_n::kappa_sample.p * (kT_loc[z].p - kT_loc[z - 1].p);
					activeField[z].p += -(qInterface.p - qBulk.p) * (deltaTime / c_ph);
				}
				else
                {
					activeField[z].p += ttm_n::kappa_sample.p * (kT_loc[z + 1].p - 2.0 * kT_loc[z].p + kT_loc[z - 1].p) * (deltaTime / c_ph);
				}	
			}
			
            // rescale substrate cells
            qInterface.e /= this->cellScale_sub;
            qInterface.p /= this->cellScale_sub;

            // substrate loop
			for (int z = this->dim + 1; z < this->dim + this->dim_sub + 1; z++)
            {

				attenuation = (ttm_n::delta_sample / ttm_n::delta_substrate) *
                    exp(-(this->dim / ttm_n::delta_sample) - ((double)(z - this->dim - 1.0) * this->cellScale_sub / ttm_n::delta_substrate));
                c_el = ttm_n::gamma_substrate * kT_loc[z].e / SI_kBoltzmann;
                c_ph = ttm_n::c_ph_substrate * ttm_n::funcEinstein(kT_loc[z].p, ttm_n::Debye_substrate);
                g_ep = ttm_n::g_ep_substrate;

				//////////////
				// ELECTRONS
				//////////////

                activeField[z].e += -g_ep* (kT_loc[z].e - kT_loc[z].p) * (deltaTime / c_el);
				activeField[z].e += ttm_n::funcLaserPower(time) * attenuation * SI_kBoltzmann * (deltaTime / c_el);
				
                kappa = ttm_n::kappa_substrate.e * (kT_loc[z].e / kT_loc[z].p);

                // special treatment for non-local terms at interface
				if (z == (this->dim + 1))
                {
					qBulk.e = -kappa * (kT_loc[z + 1].e - kT_loc[z].e);
					activeField[z].e += -(qBulk.e - qInterface.e) * (deltaTime / c_el);
				}
				else
                {
					activeField[z].e += kappa * (kT_loc[z + 1].e - 2.0 * kT_loc[z].e + kT_loc[z - 1].e) * (deltaTime / c_el);
					gradKappa = 0.5f * ttm_n::kappa_substrate.e * (kT_loc[z + 1].e / kT_loc[z + 1].p - kT_loc[z - 1].e / kT_loc[z - 1].p);
					activeField[z].e += 0.5f * gradKappa * (kT_loc[z + 1].e - kT_loc[z - 1].e) * (deltaTime / c_el);
				}
			
				//////////////
				// PHONONS
				//////////////
			
				activeField[z].p += g_ep * (kT_loc[z].e - kT_loc[z].p) * (deltaTime / c_ph);

                // special treatment for non-local terms at interface
				if (z == (this->dim + 1))
                {
					qBulk.p = -ttm_n::kappa_substrate.p * (kT_loc[z + 1].p - kT_loc[z].p);
					activeField[z].p += -(qBulk.p - qInterface.p) * (deltaTime / c_ph);
				}
				else
                {
					activeField[z].p += ttm_n::kappa_substrate.p * (kT_loc[z + 1].p - 2.0 * kT_loc[z].p + kT_loc[z - 1].p) * (deltaTime / c_ph);
				}

			}

			// at predictor step, switch fields
            // at corrector step, write result to output
			if( !integratorStep )
            {
				for(int z = 1; z < this->dim + this->dim_sub + 1; z++)
                {
					kT_loc[z].e += ttm_n::predictor[z].e;
					kT_loc[z].p += ttm_n::predictor[z].p;
				}
				activeField = ttm_n::corrector;
				this->time += deltaTime;
			} 
			else
            {
                for (int z = 1; z < this->dim + this->dim_sub + 1; z++)
                {
					kT_loc[z].e = ttm_n::result[z].e + 0.5 * (ttm_n::predictor[z].e + ttm_n::corrector[z].e);
					kT_loc[z].p = ttm_n::result[z].p + 0.5 * (ttm_n::predictor[z].p + ttm_n::corrector[z].p);
				}
			}
		}
	}
	
    return status;
}


int ttm_t::setFieldAndParameters(const int argc, const char *argv[])
{
	int status = STATUS_SUCCESS;
	int error;
	
	char message[1024];

    // list of all parameters accessible via cmd line
    double *argListVar_double[] = { 
        &ttm_n::power, 
        &ttm_n::delay,
        &ttm_n::tau,
        &ttm_n::gamma_sample, 
        &ttm_n::gamma_substrate,
        &ttm_n::c_ph_sample, 
        &ttm_n::c_ph_substrate,
        &ttm_n::g_ep_sample, 
        &ttm_n::g_ep_substrate,
        &ttm_n::delta_sample, 
        &ttm_n::delta_substrate,
        &ttm_n::kappa_sample.e,
        &ttm_n::kappa_substrate.e,
        &ttm_n::kappa_sample.p,
        &ttm_n::kappa_substrate.p,
        &ttm_n::kappa_interface.e,
        &ttm_n::kappa_interface.p,
        &ttm_n::Debye_sample,
        &ttm_n::Debye_substrate,
        nullptr
    };

    // list of all cmd line arguments
    const char *argListString_double[] = {
        "-2tm-power", 
        "-2tm-delay",
        "-2tm-tau",
        "-2tm-gamma", 
        "-2tm-substrate-gamma",
        "-2tm-c.p",
        "-2tm-substrate-c.p",
        "-2tm-gep", 
        "-2tm-substrate-gep",
        "-2tm-delta", 
        "-2tm-substrate-delta",
        "-2tm-kappa.e",
        "-2tm-substrate-kappa.e",
        "-2tm-kappa.p", 
        "-2tm-substrate-kappa.p",
        "-2tm-kapitza.e", 
        "-2tm-kapitza.p",
        "-2tm-Debye",
        "-2tm-substrate-Debye",
        "\0"
    };

    // list of each cmd line arguments's unit
    const double argListUnit_double[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        SI_milliElectronvolt, SI_milliElectronvolt, 0.0 };

    // go through cmd line args
    int argCounter = 0;
    while (argListVar_double[argCounter] != nullptr)
    {
        error = getArg(argListVar_double[argCounter], argListString_double[argCounter], argc, argv);
        if ((error & ~STATUS_SUCCESS) && (error & ~STATUS_ARG_NOT_FOUND))
        {
            sprintf(message, "temperature initialization error '%s' for '%s'", 
                getErrorString(error), 
                argListString_double[argCounter]);
            printError(message);
            abortProgram(status | error);
        }
        else if (error == STATUS_SUCCESS)
        {
            *argListVar_double[argCounter] *= argListUnit_double[argCounter];
        }
        status |= error;
        argCounter++;
    }
    
    // scale transport parameters with cell size
    ttm_n::kappa_interface.e /= this->cellLength;
    ttm_n::kappa_interface.p /= this->cellLength;
    ttm_n::kappa_sample.e /= sqr(this->cellLength);
    ttm_n::kappa_sample.p /= sqr(this->cellLength);
    ttm_n::kappa_substrate.e /= sqr(this->cellLength * this->cellScale_sub);
    ttm_n::kappa_substrate.p /= sqr(this->cellLength * this->cellScale_sub);

    ttm_n::predictor = (double2tm *)safeMalloc(this->size * sizeof(double2tm));
    ttm_n::corrector = (double2tm *)safeMalloc(this->size * sizeof(double2tm));
    ttm_n::result    = (double2tm *)safeMalloc(this->size * sizeof(double2tm));
   
    // set initial condition (homogeneous)
    for (int z = 0; z < this->size; z++)
    {
        this->kT_field_64bit[z].e = this->kT_eq;
        this->kT_field_64bit[z].p = this->kT_eq;
    }
    this->kT.e = this->kT_eq;
    this->kT.p = this->kT_eq;

	return status;
}


void ttm_t::show()
{

	printf("temperature model '%s'\n", __TTM_FUNC__);
    printf("grid parameters\n");
    printf("  sample dim %d\n", this->dim);
    printf("  substrate dim %d\n", this->dim_sub);
    printf("  substrate scale %d\n", this->cellScale_sub);
    printf("  axis '%s'\n", getAxisString(TTM_AXIS));
    printf("laser parameters\n");
    printf("  delay = %g ps\n", ttm_n::delay / SI_picoSecond);
    printf("  width = %g ps\n", ttm_n::tau / SI_picoSecond);
    printf("  power = %g W/m^3\n", ttm_n::power);
    printf("sample parameters\n");
    printf("  el. spec. heat %g J/m^3/K^2\n", ttm_n::gamma_sample);
    printf("  ph. spec. heat %g J/m^3\n", ttm_n::c_ph_sample);
    printf("  el.-ph. coupling %g W/m^3/K\n", ttm_n::g_ep_sample);
    printf("  Debye temperature %g meV\n", ttm_n::Debye_sample / SI_milliElectronvolt);
    printf("  attenuation %g layers\n", ttm_n::delta_sample);
    printf("  el. therm. cond. %g W/m/K\n", ttm_n::kappa_sample.e * sqr(this->cellLength));
    printf("  ph. therm. cond. %g W/m/K\n", ttm_n::kappa_sample.p * sqr(this->cellLength));
    printf("substrate parameters\n");
    printf("  el. spec. heat %g J/m^3/K^2\n", ttm_n::gamma_substrate);
    printf("  ph. spec. heat %g J/m^3\n", ttm_n::c_ph_substrate);
    printf("  el.-ph. coupling %g W/m^3/K\n", ttm_n::g_ep_substrate);
    printf("  Debye temperature %g meV\n", ttm_n::Debye_substrate / SI_milliElectronvolt);
    printf("  attenuation %g layers\n", ttm_n::delta_substrate);
    printf("  el. therm. cond. %g W/m/K\n", ttm_n::kappa_substrate.e * sqr(this->cellLength * this->cellScale_sub));
    printf("  ph. therm. cond. %g W/m/K\n", ttm_n::kappa_substrate.p * sqr(this->cellLength * this->cellScale_sub));
    printf("misc parameters\n");
    printf("  el. interface. cond. %g W/K/m^2\n", ttm_n::kappa_interface.e * this->cellLength);
    printf("  ph. interface. cond. %g W/K/m^2\n", ttm_n::kappa_interface.p * this->cellLength);

    return;
}


double ttm_n::funcEinstein(double kT, double Debye)
{    
	const double epsilon = 1.0e-3;

    // don't divide by zero
    if (kT / Debye < epsilon)
	{
        kT = epsilon * Debye;
    }
	
    // rough ratio between Einstein and Debye temperatures
	double beta = 0.8060 * Debye / kT;
	    
	// for small beta values return Dulong Petit
	// otherwise  return Einstein model specific heat
    if (beta < epsilon)
	{
        return 1.0;
    }
    else
	{
        return  (sqr(beta) * exp(beta) / sqr(exp(beta) - 1.0));
    }
}


double ttm_n::funcLaserPower(double time)
{
    if ( ttm_n::tau > 0.0)
    {
        return ttm_n::power * exp(-0.5 * sqr((time - ttm_n::delay) / ttm_n::delay));
    }
    else
    {
        return 0.0;
    }
}


int ttm_n::getStepSize(ttm_t ttm, double time, double deltaTime)
{

    int stepMultiple = 1;
    double2tm omega;
    double omegaMax = -1.0;

    // check sample layers for fastest dynamic
    for (int z = 1; z < ttm.dim + 1; z++)
    {
        omega.e = (ttm_n::kappa_sample.e * (ttm.kT_field_64bit[z].e / ttm.kT_field_64bit[z].p)) /
            (ttm_n::gamma_sample * ttm.kT_field_64bit[z].e / SI_kBoltzmann);
        omega.p = (ttm_n::kappa_sample.p) /
            (ttm_n::c_ph_sample * ttm_n::funcEinstein(ttm.kT_field_64bit[z].p, ttm_n::Debye_sample));

        if (omega.e > omegaMax)
        {
            omegaMax = omega.e;
        }       
        if (omega.p > omegaMax)
        {
            omegaMax = omega.p;
        }
    }

    // check substrate layers for fastest dynamics
    for (int z = ttm.dim + 1; z < ttm.dim + ttm.dim_sub + 1; z++)
    {
        omega.e = (ttm_n::kappa_sample.e * (ttm.kT_field_64bit[z].e / ttm.kT_field_64bit[z].p)) /
            (ttm_n::gamma_sample * ttm.kT_field_64bit[z].e / SI_kBoltzmann);
        omega.p = (ttm_n::kappa_sample.p) /
            (ttm_n::c_ph_sample * ttm_n::funcEinstein(ttm.kT_field_64bit[z].p, ttm_n::Debye_sample));

        if (omega.e > omegaMax)
        {
            omegaMax = omega.e;
        }
        if (omega.p > omegaMax)
        {
            omegaMax = omega.p;
        }
    }

    // reduce step size to match maximum allowed step size
    const double deltaTimeMax = 0.2 / omegaMax;
    while (deltaTime > deltaTimeMax)
    {
        deltaTime /= 2.0;
        stepMultiple *= 2;
    } 
    
    // check for zero
    if (stepMultiple > 1024)
    {
        char message[1024];
        sprintf(message, "temperature integration error, time step is suspiciously small\n");
        sprintf(message + strlen(message), "  deltaTime.2tm / deltaTime.llg = %d\n", stepMultiple);
        sprintf(message + strlen(message), "  simulation might be stuck");
        printConflict(message, __FILE__, __LINE__);
    }

    return stepMultiple;
}


#else
#error compilation failed: multiple temperature models defined!
#endif
