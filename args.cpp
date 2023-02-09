
#include "args.hpp"

void Args::readCmdArguments(int argc, char* argv[]) {

    /*
     * Function to parse the command line arguments
     * @param number of arguments 'argc', list of arguments 'argv'
     */

    int i = 1;

    //flag to set which option is already being used. to prevent using multiple times
    bool ist = false, isT = false, iso = false, isd = false, isg = false,
            ism = false, isf = false, isr = false, isp = false, isl = false,
            isnm = false, isn = false, iss = false, isS = false, isld = false,
            issw = false, isloo = false, isprec = false, isv = false,
            iswin = false, islin = false, iscont = false, istype = false,
            isbase = false, isLog = false, isor = false, isaic_parcnt=false,
            isit = false, isio = false, isid = false, ispmf = false;

    // Check the number of parameters
    if (argc < 2) {
        // Tell the user how to run the program
        cerr << Args::USAGE << endl;
        exit(-1);
    } else {
        while (i < argc) {
            string arg = argv[i];
            if ((arg == "-h") || (arg == "--help")) {
                cout << Args::USAGE << endl << VERSION << endl;
                exit(0);
            } else if ((arg == "-t") || (arg == "--train")) {
                if (ist == true) {
                    cerr << "Error:: -t --train        Train file specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                ist = true;

                g_trnFile = argv[i + 1];
                if (!Helper::is_file_exist(g_trnFile)) {
                    cerr << "Fatal Error:: No such Train file or directory: " << g_trnFile << endl;
                    cerr << "VALID OPTION: -t --train        Training Data File location. (Valid file Path)" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                i = i + 2; //go for next command line parameter
            } else if ((arg == "-T") || (arg == "--Test")) {
                if (isT == true) {
                    cerr << "Error:: -T --Test        Test file specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isT = true;

                g_tstFile = argv[i + 1];
                if (!Helper::is_file_exist(g_tstFile)) {
                    cerr << "Fatal Error:: No such Test file or directory: " << g_tstFile << endl;
                    cerr << "VALID OPTION: -T --Test        Testing Data File location. (Valid file Path)" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                i = i + 2; //go for next command line parameter
            } else if ((arg == "-log") || (arg == "--log")) {
                if (isLog == true) {
                    cerr << "Error:: -log --log        Log file specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isLog = true;
                g_logFile = argv[i + 1];
                i = i + 2; //go for next command line parameter
            } else if ((arg == "-o") || (arg == "--obj")) {
                if (iso == true) {
                    cerr << "Error:: -o --obj        Objective function specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                iso = true;
                g_objec = argv[i + 1];
                if ((g_objec == "mse") || (g_objec == "nmse") || (g_objec == "pearson") ||
                        (g_objec == "theil") || (g_objec == "med.err") ||
                        (g_objec == "med.ae") || (g_objec == "mase") || (g_objec == "ub") ||
                        (g_objec == "lb") || (g_objec == "mwe")  || (g_objec == "uwae") ||
                        (g_objec == "abserr") ||  (g_objec == "cw") ||  (g_objec == "cwy") || 
                        (g_objec == "aic")  || (g_objec == "mre") || (g_objec == "mwrae") ) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Objective Function: " << g_objec << endl;
                    cerr << "VALID OPTION: -o --obj          Objective Function (String)" << endl <<
                            "\t[Available Options: mse, nsme, pearson, theil, med.err, med.ae, mase, ub, lb, cw, cwy, aic, mwrae]" << endl <<
                            "\tDefault = mse." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-d") || (arg == "--delta")) {
                if (isd == true) {
                    cerr << "Error:: -d --delta        Delta/penalty specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isd = true;

                g_delta = atof(argv[i + 1]);
                if (Helper::is_between(g_delta, 0.0f, 1.0f) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Delta: " << g_delta << endl;
                    cerr << "VALID OPTION: -d --delta        Delta/Penalty in fitness function. (Real Number 0.0 <= d <= 1.0)" << endl <<
                            "\tDefault=0.35." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-g") || (arg == "--num-gen")) {
                if (isg == true) {
                    cerr << "Error:: -g --num-gen        Number of generations specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isg = true;
                g_numGen = atoi(argv[i + 1]);
                if (Helper::is_between(g_numGen, 1, 10000) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Generations: " << g_numGen << endl;
                    cerr << "VALID OPTION:  -g --num-gen      Number of Generations (Integer 1 <= g <= 10000)" << endl <<
                            "\tDefault = 200" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-m") || (arg == "--mut-rate")) {
                if (ism == true) {
                    cerr << "Error:: -m --mut-rate        Mutation Rate specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                ism = true;
                g_mu_rate = atof(argv[i + 1]);
                if (Helper::is_between(g_mu_rate, 0.0f, 1.0f) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Mutation Rate: " << g_mu_rate << endl;
                    cerr << "VALID OPTION:  -m --mut-rate     Mutation Rate (Real Number 0.0 <= m <= 1.0)" << endl <<
                            "\tDefault = 0.20" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-f") || (arg == "--frac-depth")) {
                if (isf == true) {
                    cerr << "Error:: -f --frac-depth        Fraction's depth specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isf = true;

                g_depth = (size_t) atoi(argv[i + 1]);
                if (Helper::is_between(g_depth, (size_t) 0, (size_t) 20) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Fraction's Depth: " << g_depth << endl;
                    cerr << "VALID OPTION:  -f --frac-depth   Depth of the Continued Fraction (Integer 0 <= f <= 20)" << endl <<
                            "\tDefault = 4" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-r") || (arg == "---reset-root")) {
                if (isr == true) {
                    cerr << "Error:: -r --reset-root        Reset Root specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isr = true;
                g_reset_stuck = atoi(argv[i + 1]);
                if (Helper::is_between(g_reset_stuck, 1, g_numGen) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Reset Root after stuck for consecutive generations: " << g_reset_stuck << endl;
                    cerr << "VALID OPTION:  -r --reset-root   Reset Root of population tree after stuck for r generations. (Integer 1 <= r <= g)" << endl <<
                            "\tDefault = 5" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-p") || (arg == "---pop-init")) {
                if (isp == true) {
                    cerr << "Error:: -p --pop-init        Population Initialisation specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isp = true;
                g_pop_init = argv[i + 1];
                string pi = argv[i + 1];
                if (pi == "avg") {
                    g_pop_init_method = 1;
                    i = i + 2; //go for next command line parameter
                } else if (pi == "sol") {
                    g_pop_init_method = 2;
                    g_solnFile = argv[i + 2];
                    if (!Helper::is_file_exist(g_solnFile)) {
                        cerr << "Fatal Error:: Population Initiliastion Method. No such solution file exists: " << g_solnFile << endl;
                        cerr << "VALID OPTION: -p --pop-init     Population Initiliastion Method." << endl <<
                                "\t[Available Options: avg (initilise linear term of solution with average of the target score)," << endl <<
                                "\t\tsol <file-name of previous solution> (Initialise population with a known solution)]" << endl <<
                                "\t\tlinreg <OPTIONAL: file-name of previous solution> (Initialise with linear regression, optionally starting from known solution)]"
                                << endl;
                        cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                        exit(-1);
                    }
                    i = i + 3; //go for next command line parameter
                } else if (pi == "linreg") {
                    g_pop_init_method = 3;
                    g_solnFile = argv[i + 2];
                    // If no solution file is provided for the linear regression option
                    if (i + 2 >= argc || g_solnFile[0] == '-') {
                        g_solnFile = (char *) "";
                        i = i + 2;
                    }// solution file is provided
                    else {
                        if (!Helper::is_file_exist(g_solnFile)) {
                            cerr << "Fatal Error:: Population Initiliastion Method. No such solution file exists: " << g_solnFile << endl;
                            cerr << "VALID OPTION: -p --pop-init     Population Initiliastion Method." << endl <<
                                    "\t[Available Options: avg (initilise linear term of solution with average of the target score)," << endl <<
                                    "\t\tsol <file-name of previous solution> (Initialise population with a known solution)]" << endl <<
                                    "\t\tlinreg <OPTIONAL: file-name of previous solution> (Initialise with linear regression, optionally starting from known solution)]"
                                    << endl;
                            cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                            exit(-1);
                        }
                        i = i + 3;
                    }
                } else {
                    cerr << "Fatal Error:: Population Initiliastion Method. No such solution file exists: " << g_solnFile << endl;
                    cerr << "VALID OPTION: -p --pop-init     Population Initiliastion Method." << endl <<
                            "\t[Available Options: avg (initilise linear term of solution with average of the target score)," << endl <<
                            "\t\tsol <file-name of previous solution> (Initialise population with a known solution)]" << endl <<
                            "\t\tlinreg <OPTIONAL: file-name of previous solution> (Initialise with linear regression, optionally starting from known solution)]"
                            << endl;
                    exit(-1);
                }
            } else if ((arg == "-l") || (arg == "--ls-gen")) {
                if (isl == true) {
                    cerr << "Error:: -l --ls-gen        Run Local Search specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isl = true;
                g_local_search = atoi(argv[i + 1]);
                if (Helper::is_between(g_local_search, 1, g_numGen) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of when to run Local Search: " << g_local_search << endl;
                    cerr << "VALID OPTION: -l --ls-gen       Run Local Search at each l(th) gen (Integer 1 <= l <= g)" << endl <<
                            "\tDefault = 1" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-lr" || arg == "--linear-prob")) {
                if (islin == true) {
                    cerr << "Error:: -lr --linear-prob        Linear regression optimization probability specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                }
                islin = true;
                g_linreg_chance = atof(argv[i + 1]);
                if (Helper::is_between(g_linreg_chance, 0.0, 1.0) == true) {
                    i = i + 2;
                } else {
                    cerr << "Fatal Error:: Invalid Value of linear regression optimization probability: " << g_linreg_chance << endl;
                    cerr << "VALID OPTION: -lr --linear-prob      Probability of running linear regression before each local optimization step." << endl <<
                            "\t(Real Number 0.0 <= lr <= 1.0) Default = 0.0" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-nm") || (arg == "--nelder-mead")) {
                if (isnm == true) {
                    cerr << "Error:: -nm --nelder-mead        Nelder-Mead Parameters specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isnm = true;
                string str(argv[i + 1]);
                string delimiter = ":";
                vector<string> v = Helper::split(str, delimiter);

                g_nm_runs = stoi(v[0]);
                g_nm_max_iters = stoi(v[1]);
                g_nm_reset_stag = stoi(v[2]);
                if ((Helper::is_between(g_nm_runs, 1, 100) == true) && (Helper::is_between(g_nm_max_iters, 1, 2500) == true) && (Helper::is_between(g_nm_reset_stag, 1, 100) == true)) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Nelder-Mead Parameters: " << str << endl;
                    cerr << "Options parsed as: nr=" << g_nm_runs << "\titer=" << g_nm_max_iters << "\treset=" << g_nm_reset_stag << endl;
                    cerr << "VALID OPTION: -nm --nelder-mead  Parameters for Nelder-Mead (NM) (all Integers): <nr>:<gen>:<reset>" << endl <<
                            "\tnr=number of runs of NM (1 to 100), gen=max generations per run of NM (1 to 1000)" << endl <<
                            "\treset=NM reset after stagnate for consecutive generations (1 to 100)" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-n") || (arg == "--norm")) {
                if (isn == true) {
                    cerr << "Error:: -n --norm        Normalisation specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isn = true;
                g_norm = atoi(argv[i + 1]);
                if (Helper::is_between(g_norm, 0, 3) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Normalisation: " << g_norm << endl;
                    cerr << "VALID OPTION: -n --norm         Data Normalisation mathod." << endl <<
                            "\t[Available Options: 0=No Normalisation, 1=MinMax, 2=z-Score, 3=Decimal Scaling.]" << endl <<
                            "\tDefault = 0" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-s") || (arg == "--seed")) {
                if (iss == true) {
                    cerr << "Error:: -s --seed        Seed for Random Number specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                iss = true;
                g_seed = (uint_fast32_t) strtoul(argv[i + 1], nullptr, 10);
                if (g_seed > 0) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value for Seed: " << g_seed << endl;
                    cerr << "VALID OPTION: -s --seed          Seed for Random Number Generator. (Positive Integer)" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                //                g_rreal = randreal(g_seed);
                //                g_rint = randint(g_seed);
            } else if ((arg == "-S") || (arg == "--is-serial")) {
                if (isS == true) {
                    cerr << "Error:: -S --is-serial        Serial Mode of Optimisation specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isS = true;
                string argv_serial = argv[i + 1];
                if (argv_serial == "true") g_run_serial = true;
                else if (argv_serial == "false") g_run_serial = false;
                else {
                    cerr << "Fatal Error:: Invalid Value for Serial Mode of Optimisation: " << argv_serial << endl;
                    cerr << "VALID OPTION: -S --is-serial         Run Local Search Optimisation in Serial Mode? (String)" << endl <<
                            "\t[Available Options: true, false.]" << endl <<
                            "\tDefault = false" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                i = i + 2; //go for next command line parameter
            } else if ((arg == "-v") || (arg == "--verbose")) {
                if (isv == true) {
                    cerr << "Error:: -v --verbose        Verbosity mode of output specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isv = true;
                g_verbose_mode = atoi(argv[i + 1]);
                if (Helper::is_between(g_verbose_mode, 0, 4) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Verbose Mode: " << g_verbose_mode << endl;
                    cerr << "VALID OPTION: -v --verbose          Verbosity mode of output. (Integer: 0 <= v <=4)" << endl <<
                            "\tDefault = 0" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if (strcmp(argv[i], "-c") == 0) {
                g_numFrac = argv[i + 1];
            } else if ((arg == "-prec") || (arg == "--out-prec")) {
                if (isprec == true) {
                    cerr << "Error:: -prec --out-prec        Set Precision in Output specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isprec = true;
                g_PREC = atoi(argv[i + 1]);
                if (Helper::is_between(g_PREC, 1, 18) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Output Precision: " << g_PREC << endl;
                    cerr << "VALID OPTION: -op --out-prec          Set Precision in Output. (1 <= prec <= 18)" << endl <<
                            "\tDefault = 6" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-loo") || (arg == "--loo-cv")) {
                if (isloo == true) {
                    cerr << "Error:: -loo --loo-cv        Leave-one-out Cross-validation specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isloo = true;
                g_loocv = true;
                i = i + 1;
            } else if ((arg == "-id") || (arg == "--inc-depth")) {
                if (isid == true) {
                    cerr << "Error:: -id --inc-depth      Incremental Depth switch specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isid = true;
                g_incremental_depth = true;
                i = i + 1; //go for next command line parameter
            } 
            // -pmf --pr-model-mf Enabling the usage of model form previous Depth
            else if ((arg == "-pmf") || (arg == "----pr-model-mf")) {
                if (ispmf == true) {
                    cerr << "Error:: -pmf --pr-model-mf      Usage of Previous Models as Meta-Feature switch specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                ispmf = true;
                g_prev_model_as_var = true;
                i = i + 1; //go for next command line parameter
            } 
            else if ((arg == "-ld") || (arg == "--ls-data")) {
                if (isld == true) {
                    cerr << "Error:: -ld --ls-data      Percentage of Samples used to eval model in Local Search specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isld = true;
                g_percnt_data = atoi(argv[i + 1]);
                if (Helper::is_between(g_percnt_data, 1, 100) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Percentage of Data for Evaluation of LS: " << g_percnt_data << endl;
                    cerr << "VALID OPTION: -ld --ls-data          Percentage of Samples used to eval model in Local Search. (1 <= ld <= 100)" << endl <<
                            "\tDefault = 100" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-loo") || (arg == "--loo-cv")) {
                if (isloo == true) {
                    cerr << "Error:: -loo --loo-cv        Leave-one-out Cross-validation specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isloo = true;
                g_loocv = true;
                i = i + 1;
            }
            else if (arg == "--cont") {
                if (iscont == true) {
                    cerr << "Error:: --cont        Continuous Eval specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                iscont = true;
                string argv_cont = argv[i + 1];
                if (argv_cont == "true") g_cont = true;
                else if (argv_cont == "false") g_cont = false;
                else {
                    cerr << "Fatal Error:: Invalid Value for Continuous Eval: " << argv_cont << endl;
                    cerr << "VALID OPTION: --cont         Continuous Eval (String)" << endl <<
                            "\t[Available Options: true, false.]" << endl <<
                            "\tDefault = false" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                i = i + 2; //go for next command line parameter
            } else if ((arg == "-w") || (arg == "--window")) {
                if (iswin == true) {
                    cerr << "Error:: -w --window      Set the size of window (in percentage of training data) to use in local search specified multiple times ." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                iswin = true;
                int win = atoi(argv[i + 1]);
                if (Helper::is_between(win, 1, 100) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Percentage of Data Window for Evaluation of LS: " << g_window_size << endl;
                    cerr << "VALID OPTION: -w --window      Set the size of window (in percentage of training data) to use in local search." << endl <<
                            "\tDefault = 100" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                g_window_size = win / 100.0; //convert into double value of the percentage
            } else if ((arg == "-y") || (arg == "--type")) {
                if (istype == true) {
                    cerr << "Error:: -y --type        Type of output specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                istype = true;
                g_type = atoi(argv[i + 1]);
                if (Helper::is_between(g_type, 0, 3) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Value of Type: " << g_type << endl;
                    cerr << "VALID OPTION: -y --type          Type of output. (Integer: 0 <= y <= 3)" << endl <<
                            "\tDefault = 0" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-b") || (arg == "--base")) {
                if (isbase == true) {
                    cerr << "Error:: -b --base        Base specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isbase = true;
                g_base = atof(argv[i + 1]);
                if (Helper::is_between(g_base, 0, 10000) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Base: " << g_base << endl;
                    cerr << "VALID OPTION: -b --base          Type of output. (Integer: 0 <= b <= 10000)" << endl <<
                            "\tDefault = 0" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            } else if ((arg == "-sw") || (arg == "--samp-weight")) {
                if (issw == true) {
                    cerr << "Error:: -sw --samp-weight      Sample weight flag specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                issw = true;
                g_weight_present = true;
                i += 1;

            } else if ((arg == "-or") || (arg == "--origin")) {
                if (isor == true) {
                    cerr << "Error:: -or --origin      Origin flag specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isor = true;
                g_origin = true;
                i += 1;

            } else if ((arg == "-pc") || (arg == "--param-count")) {
                if (isaic_parcnt == true) {
                    cerr << "Error:: -pc --param-count     AIC meta-feature parameter counter specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isaic_parcnt = true;
                g_aic_paramCount = true;
                g_paramCountFile = argv[i + 1];
                i += 2;

            } else if ((arg == "-io") || (arg == "--integer-only")) {
                if (isio == true) {
                    cerr << "Error:: -io --integer-only   Fag specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isio = true;
                g_int_only = true;
                i += 1;

            } else if ((arg == "-it") || (arg == "--integer-terms")) {
                if (isit == true) {
                    cerr << "Error:: -it --integer-terms      Flag specified multiple times." << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
                isit= true;
                g_int_terms = atoi(argv[i + 1]);
                if (Helper::is_between(g_int_terms, 0, 3) == true) {
                    i = i + 2; //go for next command line parameter
                } else {
                    cerr << "Fatal Error:: Invalid Flag: " << g_base << endl;
                    cerr << "VALID OPTION: -it --integer-terms          Type of output. (Integer: 0 <= it <= 3)" << endl <<
                            "\tDefault = 0" << endl;
                    cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
                    exit(-1);
                }
            }  else if ((arg == "-V") || (arg == "--version")) {
                cout << "Continued Fraction Regression." << endl << VERSION << endl;
                exit(0);
            } else {
                i = i + 1;
            }
        }

        //check for mwe objective function, if -sw was supplied
        if ((g_objec == "mwe" ) && issw != true) {
            cerr << "Fatal Error:: Sample weight (-sw or -sw --samp-weight) is required for Mean Weighetd Error (mwe) objective fucntions." << endl;
            exit(-1);
        }

        // Check test file not provided with LOOCV
        if ( g_loocv && strlen(g_tstFile) > 0) {
            cerr << "Fatal Error:: Testing Data is not allowed in Leave-One-Out Cross-Validation." << endl;
            cerr << "You can run \'./bin/main --help\' to view all available options." << endl;
            exit(-1);
        }

    }

    // Generate random seed if not provided one
    if(g_seed==0){
        randint s_rint = randint();
        g_seed = s_rint(1, numeric_limits<int>::max());
    }

    // Debug arguments
    if( g_verbose_mode > 0) {

        // Debug Args Vs Parameters
        cout << "Cmd: " << argv[0];
        for (int pi = 1; pi < argc; pi++) {

            // Start next argument on a new line
            if( argv[pi][0] == '-' )    
                cout << endl << "\t";
            
            // Write flag or argument
            cout << argv[pi] << " ";
        }
        cout << endl;
    }
}

void Args::printProgramParameters() {
    /*
     * Print Parameters used in the program
     */

    cout << endl;
    cout << "======================== PARAMETER VALUE =====================" << endl;
    cout << "-t --train         " << g_trnFile << endl;
    cout << "-T --Test          " << g_tstFile << endl;
    cout << "-log               " << g_logFile << endl;
    cout << "-o --obj           " << g_objec << endl;
    if (g_aic_paramCount==true)       cout << "-pc --param-count " << g_paramCountFile << endl;
    cout << "-d --delta         " << g_delta << endl;
    cout << "-g --num-gen       " << g_numGen << endl;
    cout << "-m --mut-rate      " << g_mu_rate << endl;
    cout << "-f --frac-depth    " << g_depth << endl;
    cout << "-r --reset-root    " << g_reset_stuck << endl;
    if (g_pop_init_method == 0)         cout << "-p --pop-init      " << "Random" << endl;
    else if (g_pop_init_method == 1)    cout << "-p --pop-init      " << "avg" << endl;
    else if (g_pop_init_method == 2)    cout << "-p --pop-init      " << "sol " << g_solnFile << endl;
    else if (g_pop_init_method == 3)    cout << "-p --pop-init      " << "linreg " << g_solnFile << endl;
    cout << "-l --ls-gen        " << g_local_search << endl;
    cout << "-lr --linear-prob  " << g_linreg_chance << endl;
    cout << "-nm --nelder-mead  " << g_nm_runs << ":" << g_nm_max_iters << ":" << g_nm_reset_stag << endl;
    cout << "-n --norm          " << g_norm << endl;
    cout << "-S --is-serial     " << boolalpha << g_run_serial << endl;
    if (g_seed > 0) cout << "-s --seed          " << g_seed << endl;
    else cout << "-s --seed          " << "random_device" << endl;
    cout << "-ld --ls-data      " << g_percnt_data << endl;
    if (g_prev_model_as_var == true)    cout << "-pmf --pr-model-mf " << g_prev_model_as_var << endl;
    if (g_incremental_depth == true)    cout << "-id --inc-depth    " << g_incremental_depth << endl;
    cout << "-c --cont          " << boolalpha << g_cont <<endl;
    cout << "-y --type          " << g_type << endl;
    cout << "-b --base          " << g_base << endl;
    cout << "-w --window        " << g_window_size << endl;
    if (g_weight_present==true)       cout << "-sw --samp-weight  present" << endl;
    else                            cout << "-sw --samp-weight  not present" << endl;
    cout << "-loo --loo-cv      " << boolalpha << g_loocv << endl;
    cout << "-prec --out-prec   " << g_PREC << endl;
    cout << "-v --verbose       " << g_verbose_mode << endl;
    cout << "-or --origin       " << g_origin << endl;
    cout << "============================================================" << endl;
}

const string Args::USAGE =

    R"(Continued Fraction Regression.

    REQUIRED ARGUMENTS:

        -t --train         Training Data File location.

    OPTIONS:

        -b --base           If the representation for the model is CCL, selects the base used (Integer 0 < b <= 10000)
                            Default = 2

        --cont              Data used to eval model in Local Search is continuous instead of randomly selected (String)
                            [Available Options: true, false]
                            Default = false

        -d --delta          Delta/Penalty in fitness function. (Real Number 0.0 <= d <= 1.0)
                            Default=0.35.

        -f --frac-depth     Depth of the Continued Fraction (Integer 0 <= f <= 20)
                            Default = 4

        -g --num-gen        Number of Generations (Integer 1 <= g <= 10000)
                            Default = 200

        -h --help           Show this screen.

        -io,                Only allow integer coefficients (no argument)
        --interger-only

        -it,                Force the terms f(x) = {gi(x), hi(x)} to resolve to integer values. An integer argument must be provided:
        --integer-terms         0 = integers are not forced
                                1 = floor(f(x)) is take
                                2 = round(f(x)) is taken
                                3 = ceil(f(x)) is taken

                            Local search behaviour may be modified based on this parameter to focus on the exploration of the natural numbers 
                            instead of the reals

        -l --ls-gen         Run Local Search at each l(th) gen (Integer 1 <= l <= g)
                            Default = 1

        -id --inc-depth     Flag (no argument) to enable Incremental Depth Approach where the model from previous depth 
                            will be used to Initialise the Population of current depth.

        -pmf --pr-model-mf  Flag (no argument) to enable the usage of model form previous Depth as a meta-feature in current depth.
                            The `-id` or `--inc-depth` must be enabled to use this feature.

        -ld --ls-data       Percentage of Samples used to eval model in Local Search. (1 <= ld <= 100)
                            Default = 100.

        -loo --loo-cv       Enable Leave-One-Out Cross Validation on the Training Data. No Parameter Values if required. By default it is disabled.

        -m --mut-rate       Mutation Rate (Real Number 0.0 <= m <= 1.0)
                            Default = 0.20

        -n --norm           Data Normalisation Method (Integer 0 <= n <=3)
                            [Available Options: 1=MinMax, 2=z-Score, 3=Decimal Scaling]
                            Default = 0 (No Normalisation).
    
        -nm --nelder-mead   Parameters for Nelder-Mead (NM) (all Integers): <nr>:<gen>:<reset>
                            nr=number of runs of NM (1 to 100), gen=max generations per run of NM (1 to 1000)
                            reset=NM reset after stagnate for consecutive generations (1 to 100).

        -o --obj            Objective Function (String)
                            [Available Options: 
                                mse, 
                                nsme, 
                                pearson, 
                                theil, 
                                med.err, 
                                med.ae, 
                                mase, 
                                ub lb, 
                                mwe, 
                                uwae,
                                cw,
                                mre,
                                aic,
                                mwrae
                                ]
                            Default = mse.

        -or --origin       Force the solution to pass through the origin point for x1,x2,..,xN = 0


        -p --pop-init       Population Initiliastion Method.
                            [Available Options: avg (Initialise linear term of solution with average of the target score),
                                                sol <file-name of previous solution> (Initialise population with a known solution)
                                                linreg <OPTIONAL: file-name of previous solution> (Initialise with linear regression, optionally starting from known solution)]

        -prec --out-prec    Set Precision in Output. (1 <= op <= 18)
                            Default = 6.

        -r --reset-root     Reset Root of population tree after stuck for r generations. (Integer 1 <= r <= g)
                            Default =  5

        -s --seed           Seed for Random Number Generator. (Positive Integer)
                            Default = 13.

        -sw --samp-weight   Flag (no argument) to indicate the last column is sample weight
                            The variable count in the dataset file must be (feature count)+1 when weight supplied e.g.

                            Three features      Three features + sample weight 
                            y,f1,f2,f3          y,f1,f2,f3,w
                            1,0,0,1             1,0,0,1,5
                            2,0,1,0             2,0,1,0,3
                            3,0,1,1             3,0,1,1,5

                            NB  The test data has no purpose for weights, but test and training files are handled the same 
                            Currently, you must include the weight column ```w``` in the test .csv file.
                            The code will ignore this column for the testing data

        -pc --param-count  AIC: path to the file containing the parameter-counter when meta-feature is present.
        -loo --loo-cv      Enable Leave-One-Out Cross Validation on the Training Data. No Parameter Values if required. By default it is disabled.

        -S --is-serial      Run Local Search Optimisation in Serial Mode? (String)
                            [Available Options: true, false]
                            Default = true.
        
        -t --train          Training Data File location. (Valid file Path)

        -T --Test           Testing Data File location. (Valid file Path)
        
        -v --verbose        Verbosity mode of output. (Integer: 0 <= v <=4)
                            Default = 0

        -V --version        Show version.

        -w --window         Set the size of window (in percentage of training data) to use in local search. (1 <= w <= 100)
                            Default = 100.

        -y --type           Type of representation for model output (Integer)
                            [Available Options: 0 for Continued Fraction (CF), 1 for Coefficient Continued Logarithm (CCL)]
                            Default = 0
        
    USAGE:
        .bin/main -t <train-file>
        .bin/main -t <train-file> -T <test-file>
        .bin/main -t <train-file> -T <test-file> -s 107
        .bin/main -t <train-file> -T <test-file> -g 200 -m 0.1 -nm 4:200:5 -v 2
        .bin/main -t <train-file> -T <test-file> -g 200 -m 0.5 -d 0.2 -nm 4:250:10 -f 0 -p avg -v 2
        .bin/main -t <train-file> -T <test-file> -g 200 -m 0.1 -d 0.1 -nm 4:250:10 -f 4 -p sol CFR-v2-Depth-3.sol.txt -v 2
        .bin/main --version
        
    )";
