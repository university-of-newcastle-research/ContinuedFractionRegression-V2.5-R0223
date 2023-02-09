
#include "pop.hpp"

// Constructors

population::population(data_store train, data_store test) :
    train_data(train), test_data(test),
    train_e(train_data), test_e(test_data) {
    /*
     * Initialise population with random value
     * @param train data, test data
     */

    num_var = train_data.num_var;
    names = train_data.names;

    root = new agent(3);
    initialize(root);
}

population::population(data_store train, data_store test, string soln_file) :
    train_data(train), test_data(test),
    train_e(train_data), test_e(test_data) {

    solfile = soln_file;

    num_var = train_data.num_var;
    names = train_data.names;

    root = new agent(3);
    initialize(root, soln_file);
}


population::~population() {
    /*
     * Destructor of Population Object
     */ 
    delete root;
}

void population::breed(agent* leader) {
    /*
     * Recombination Function for the Memetic Algorithm
     * @param a pointer to the agent to apply the operator
     */

    mutate(leader);

    // reached a leaf
    if (leader->depth <= 1) return;

    size_t n = leader->degree;
    recombine(leader, leader->children[n - 1]);
    recombine(leader->children[n - 1], leader);

    for (size_t i = 0; i < n - 1; ++i) {
        recombine(leader->children[i], leader->children[i + 1]);
    }

    for (size_t i = 0; i < n; ++i) {
        breed(leader->children[i]);
    }

    if (!check(root)) {
        cerr << "Recombination!\n";
        assert(false);
    }
}

string createOutFile(string fileName) {
    /*
     * Create name of the output file by appending "out.csv" at the end of file name
     * @param file name
     * @return name of the output file
     */
    string outFile = fileName + string(".D-") + to_string(g_depth) + string(".out.csv");
    return outFile;
}

void population::run(MatrixContinuedFraction& sol) {
    /*
     * Run the Memetic Algorithm Framework for the Continued Fraction Regression
     */
    double best = root->fitness[0];
    double cur_best = root->fitness[0];

    sol = root->member[0];
    
    // Creating the output file names form training & testing file names
    string buftrn(g_trnFile);
    buftrn = createOutFile(buftrn);

    string buftst(g_tstFile);
    buftst = createOutFile(buftst);

    cout << "Starting the Memetic Algorithm" << endl;

    // set output precision
    cout.precision(g_PREC);
    auto t_start = chrono::system_clock::now();

    int staleCount = 0;
    // int switch_obj_now = (g_numGen/2)+1;

    for (int i = 1; i <= g_numGen; ++i) {

        breed(root);        // Mutate and recombine the agent solutions

        // !!! Needs justification
        if (i % g_local_search == 0) {
            // diversify(root);
            local_search(root);
            if (!check(root)) {
                cerr << i << " Local search!\n";
                print_population();
                assert(false);
            }
        }

        // count the num of generations for which the fitness did not improved
        if (root->fitness[0] >= cur_best) ++staleCount;
        else staleCount = 0;

        // Reset best solution if it has not improved in g_reset_stuck iterations
        
        if (staleCount > g_reset_stuck) {

            // Randomly initialise a new solution
            if( solfile != "") {
                root->member[0] = MatrixContinuedFraction(solfile);
            } else {
                root->member[0] = MatrixContinuedFraction(g_depth, num_var);
            }

            
            eval_fit(root, 0);
            root->movedown_pocket();
            staleCount = 0;

            // !!! Needs justification
            if (!check(root)) {
                cerr << "Old-age!\n";
                print_population();
                assert(false);
            }
        }
       
        // Determine if the new generation has a better solution
        cur_best = root->fitness[0];
        if (cur_best < best) {

            // Update the best solution
            best = cur_best;
            sol = root->member[0];

            // Log each new solition to file
            dbg(fileModels,to_string(i)+",",0,0);
            if(g_verbose_mode > 0) sol.write_cfr(cout, names, Numpy, true);     // Moh: Keep vebose 0 as minimal, show best solution for v>0
            if(g_verbose_mode > 1)  sol.write_term_table(cout, names, Numpy);    // Moh: I think it should be introduced to more higher level of verbosity than 1
            dbgnl(fileModels);
        }

        // print the current generation
        if(g_verbose_mode == 0)
            cout << setw(4) << i << "\t" << setw(g_PREC+3) << sol.fitness << endl;
        else
            cout << setw(4) << i << "\t fitness:" << setw(g_PREC+3) << sol.fitness << endl;
    }
    cout << endl;

    auto t_end = chrono::system_clock::now();
    chrono::duration<double> t_diff = t_end - t_start;

    cout.precision(g_PREC);

    cout << "Best Solution" << endl;
    cout << "=============" << endl;
    cout << "Latex" << endl;
    sol.write_cfr(cout, names, Latex, false);
    cout << endl;
    cout << "Excel" << endl;
    sol.write_cfr(cout, names, Excel, false);
    cout << endl;
    cout << "Numpy" << endl;
    sol.write_cfr(cout, names, Numpy, false);
    cout << endl;
    int param_count = sol.sum_param_count(train_data.aic_param_count);
    cout << "# Free Parameters: " << param_count << endl;
    cout << endl;
    // Show the rational number representation of the coefficients and constants
    cout << "============= Rational Representation of Coefficients and Constant =============" << endl;
    cout << "Latex" << endl;
    sol.write_cfr(cout, names, Latex, false, true);
    cout << endl;
    cout << "Excel" << endl;
    sol.write_cfr(cout, names, Excel, false, true);
    cout << endl;
    cout << "Numpy" << endl;
    sol.write_cfr(cout, names, Numpy, false, true);
    cout << endl << endl;

    cout << "==================================================" << endl << endl;
    cout << "BestFitness\t " << best << "\tSolution_Size\t" << sol.used_features() << "\tDelta_Penalty\t" << sol.penalty << "\tDepth\t" << sol.depth<< endl << endl;
    print_score_header();
    cout << "Train\t";
    train_e.compute_Scores(sol);
    cout << "Test\t";
    test_e.compute_Scores(sol);


    string logFile  = g_outPath + "Run"+to_string(g_run_no-1)+".log";
    ofstream log(logFile);
    if (log.is_open()) {
        
        vector<string> xNames = vector<string>();
        for( size_t i = 0; i < names.size(); i++) {
            xNames.push_back("x"+to_string(i));
        }

        log << "Seed,Train MSE,Test MSE,Dur,Model,Train X2,Train X2 %" << endl;
        log << g_seed << "," << train_e.compute_MSE(sol) << "," << test_e.compute_MSE(sol) << "," <<  t_diff.count() << ",";
        sol.write_cfr(log, xNames, Numpy, false);
        if( g_objec == "cw") {
            log <<  "," << train_e.compute_ChiSquared(sol) << "," << train_e.compute_ChiSquaredPct(sol);
            log <<  "," << test_e.compute_ChiSquared(sol) << "," << test_e.compute_ChiSquaredPct(sol);
        }
        
        log << endl;
        log.close();
    }

    //write the solution to file ( currently ignoring the sOut/soutFile, the file name is being created inside the export_cfr fucntion)
    string soutFile = g_logFile;
    const char *sOut = soutFile.c_str();
    sol.export_cfr(sOut, g_run_no-1, names);
    train_data.write_csv_model("", g_run_no-1, true, sol);
    train_data.write_csv_model("", g_run_no-1, false, sol);

    //write train and test predictions
    if (buftrn.compare(buftst) != 0) { //check if different train and test
        const char *trnOut = buftrn.c_str();
        const char *tstOut = buftst.c_str();

        train_e.print_val(trnOut, sol);
        test_e.print_val(tstOut, sol);
    } else { // if same file is used as train and test
        const char *trnOut = buftrn.c_str();
        train_e.print_val(trnOut, sol);
    }
    if (g_norm > 0) { //write normalised test file
        string buftstNorm(g_tstFile);
        string outFileNorm = buftstNorm + ".D-" + to_string(sol.depth) + ".out.norm.csv";
        const char *tstOutNorm = outFileNorm.c_str();

        test_e.print_val_nd(tstOutNorm, sol, train_data);
    }

    g_fitness = sol.error;
}

inline void population::print_score_header() {   

    /*
     * Print the header of the scores.
     *  The measure used in the fitness function is highlighted with *
     */

    // size_t pad = g_PREC+5;

    // Column one contains "Train" or "Test" so we offset buy this amount
    // cout << "\t";

    cout << "\t";
    if (g_objec == "mse")       cout << "*mse";
    else                        cout << "mse";

    cout << "\t";
    if (g_objec == "nmse")      cout << "*nmse";
    else                        cout << "nmse";

    cout << "\t";
    if (g_objec == "pearson")   cout << "*pearson";
    else                        cout << "pearson";

    cout << "\t";
    if (g_objec == "theil")     cout << "*theil";
    else                        cout << "theil";

    cout << "\t";
    if (g_objec == "med.err")   cout << "*med.err";
    else                        cout << "med.err";

    cout << "\t";
    if (g_objec == "med.ae")    cout << "*med.ae";
    else                        cout << "med.ae";

    cout << "\t";
    if (g_objec == "mase")      cout << "*mase";
    else                        cout << "mase";

    cout << "\t";
    if (g_objec == "cw")        cout << "*Chi^2";
    else                        cout << "Chi^2";

    cout << "\t" << "% in dy range" << endl;

}

pair <double, double> population::run(int k) {
 
    /*
     * Run the Memetic Algorithm Framework for the Continued Fraction Regression for
     */

    double best = root->fitness[0];
    MatrixContinuedFraction sol = root->member[0];
    double cur_best = root->fitness[0];

    int counter = 0;

    for (int i = 1; i <= g_numGen; ++i) {

        breed(root); // doing mutation and reco
        if (i % g_local_search == 0) {
            local_search(root);

            if (!check(root)) {
                cerr << i << " Local search!\n";
                print_population();
                assert(false);
            }
        }

        if (root->fitness[0] >= cur_best) ++counter;
        else counter = 0;

        // kill if the best solution got stuck
        if (counter > g_reset_stuck) {
            //moh:mat
            root->member[0] = MatrixContinuedFraction(g_depth, num_var); //fraction(num_var);
            eval_fit(root, 0);
            root->movedown_pocket();
            counter = 0;

            if (!check(root)) {
                cerr << "Old-age!\n";
                print_population();
                assert(false);
            }
        }
        //moh: printing the entire population
        if (g_verbose_mode == 4) {
            print_population();
        }

        cur_best = root->fitness[0];
        if (cur_best < best) {
            best = cur_best;
            sol = root->member[0];
        }
    }
    
    
    pair <double, double> res_mse;
    double loocv_mse = test_e.compute_MSE(sol);
    res_mse.first = best;
    res_mse.second = loocv_mse;
    return res_mse;
}

bool population::check(agent* a) {
    if (a->fitness[0] > a->fitness[1]) {
        cerr << "Pocket! " << a->depth << ' ' << a->fitness[0] << ' ' << a->fitness[1] << endl;
        return false;
    }

    if (a->depth <= 1) return true;

    for (size_t i = 0; i < a->degree; ++i) {
        if (!check(a->children[i])) return false;

        if (a->fitness[0] > a->children[i]->fitness[0]) {
            cerr << "Parent! " << a->depth << ' ' << a->fitness[0] << ' ' << a->children[i]->fitness[0] << endl;
            return false;
        }
    }

    return true;
}

void population::print_population_coeff() {
    /*
     * To Print the coefficients of agents in the population 
     */

    deque<agent*> q;
    q.push_back(root);

    while (!q.empty()) {
        agent* a = q.front();
        q.pop_front();

        if (a->depth > 1)
            for (size_t i = 0; i < a->degree; ++i)
                q.push_back(a->children[i]);

        if (q.empty() || q.front()->depth < a->depth)
            cout << '\n';
    }
}

string population::print_agent_vars(agent* a, int sol) {
    /*
     * To Print the agent's variable selection information
     * @param pointer to the agent, pocket/current solution
     * @return string for printing
     */

    char buffer[100];
    char type[3] = "PC"; // 0=P and 1=C

    sprintf(buffer, "[%c: %6.3f: ", type[sol], a->fitness[sol]);
    string out(buffer);

    string vars = "";
    for (unsigned int i = 0; i < a->member[sol].vars; i++) {
        vars += to_string(a->member[sol].feature_active[i]);
    }
    vars += "]";
    return (out + vars);
}

void population::print_population() {
    /*
     * TO Print the entire population
     */

    deque<agent*> q;
    q.push_back(root);

    while (!q.empty()) {
        agent* a = q.front();
        q.pop_front();

        string out = print_agent_vars(a, 0);
        out += "\n";
        out += print_agent_vars(a, 1);
        cout << out << endl;

        if (a->depth > 1)
            for (size_t i = 0; i < a->degree; ++i)
                q.push_back(a->children[i]);

        if (q.empty() || q.front()->depth < a->depth)
            cout << '\n';
    }
}
