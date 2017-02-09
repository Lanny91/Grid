/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid

 Source file: tests/hadrons/Test_hadrons_rarekaon.cc

 Copyright (C) 2017

 Author: Andrew Lawson <andrew.lawson1991@gmail.com>

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

 See the full license in the file "LICENSE" in the top level distribution
 directory.
 *******************************************************************************/

#include "Test_hadrons.hpp"

using namespace Grid;
using namespace Hadrons;


int main(int argc, char *argv[])
{
    // parse command line //////////////////////////////////////////////////////
    std::string configStem;
    
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <configuration filestem> [Grid options]";
        std::cerr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    configStem = argv[1];
    
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    unsigned int light   = 0;
    unsigned int strange = 1;
    unsigned int charm   = 2;
    std::vector<double>       mass    = {.01, .04, .2};
    std::vector<std::string>  flavour = {"l", "s", "c"};
    std::vector<std::string>  solvers = {"CG_l", "CG_s", "CG_c"};
    std::string               kmom    = "0. 0. 0. 0.";
    std::string               pmom    = "1. 0. 0. 0.";
    std::string               qmom    = "-1. 0. 0. 0.";
    std::string               mqmom   = "1. 0. 0. 0.";
    std::vector<unsigned int> tKs     = {0, 16};
    unsigned int              dt_pi   = 16;
    std::vector<unsigned int> tJs     = {8};
    unsigned int              n_noise = 0;
    unsigned int              nt      = 32;

    // Global parameters.
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 1500;
    globalPar.trajCounter.end      = 1520;
    globalPar.trajCounter.step     = 20;
    globalPar.seed                 = "1 2 3 4";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // gauge field
    MGauge::Load::Par gaugePar;
    gaugePar.file = configStem;
    application.createModule<MGauge::Load>("gauge", gaugePar);
    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 12;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);

        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action   = "DWF_" + flavour[i];
        solverPar.residual = 1.0e-8;
        application.createModule<MSolver::RBPrecCG>(solvers[i],
                                                    solverPar);
    }

    // Create noise propagators for loops.
    if (n_noise > 0)
    {
        std::vector<std::string> noiseProps_l(n_noise);
        std::vector<std::string> noiseProps_c(n_noise);
        MSource::Z2::Par noisePar;
        noisePar.tA = 0;
        noisePar.tB = nt - 1;
        for (unsigned int nn = 0; nn < n_noise; ++nn)
        {
            std::string eta      = INIT_INDEX("noise", nn);
            std::string loop_l   = INIT_INDEX("loop_l", nn);
            std::string loop_c   = INIT_INDEX("loop_c", nn);

            application.createModule<MSource::Z2>(eta, noisePar);
            makePropagator(application, loop_l, eta, solvers[light]);
            makePropagator(application, loop_c, eta, solvers[charm]);
            noiseProps_l.push_back(loop_l);
            noiseProps_c.push_back(loop_c);
        }
    }

    // Translate rare kaon decay across specified timeslices.
    for (unsigned int i = 0; i < tKs.size(); ++i)
    {
        // Zero-momentum wall source propagators for kaon and pion.
        unsigned int tK     = tKs[i];
        unsigned int tpi    = (tK + dt_pi) % nt;
        std::string q_Kl_0  = INIT_INDEX("Q_l_0", tK);
        std::string q_pil_0 = INIT_INDEX("Q_l_0", tpi);
        MAKE_WALL_PROP(tK, q_Kl_0, solvers[light]);
        MAKE_WALL_PROP(tpi, q_pil_0, solvers[light]);

        // Wall sources for kaon and pion with momentum insertion. If either
        // p or k are zero, or p = k, re-use the existing name to avoid 
        // duplicating a propagator.
        std::string q_Ks_k  = INIT_INDEX("Q_Ks_k", tK);
        std::string q_Ks_p  = INIT_INDEX((kmom == pmom) ? "Q_Ks_k" : "Q_Ks_p", tK);
        std::string q_pil_k = INIT_INDEX((kmom == ZERO_MOM) ? "Q_l_0" : "Q_l_k", tpi);
        std::string q_pil_p = INIT_INDEX((pmom == kmom) ? q_pil_k : ((pmom == ZERO_MOM) ? "Q_l_0" : "Q_l_p"), tpi);
        MAKE_3MOM_WALL_PROP(tK, kmom, q_Ks_k, solvers[strange]);
        MAKE_3MOM_WALL_PROP(tK, pmom, q_Ks_p, solvers[strange]);
        MAKE_3MOM_WALL_PROP(tpi, kmom, q_pil_k, solvers[light]);
        MAKE_3MOM_WALL_PROP(tpi, pmom, q_pil_p, solvers[light]);

        // Wall sink-smeared propagators.
        /*
        std::string q_Kl_0_W  = q_Kl_0 + "_W";
        std::string q_Ks_k_W  = q_Ks_k + "_W";
        std::string q_Ks_p_W  = q_Ks_p + "_W";
        std::string q_pil_0_W = q_pil_0 + "_W";
        std::string q_pil_k_W = q_pil_k + "_W";
        std::string q_pil_p_W = q_pil_p + "_W";
        makeWallSink(application, q_Kl_0, q_Kl_0_W);
        makeWallSink(application, q_Ks_k, q_Ks_k_W, kmom);
        makeWallSink(application, q_Ks_p, q_Ks_p_W, pmom);
        makeWallSink(application, q_pil_0, q_pil_0_W);
        makeWallSink(application, q_pil_k, q_pil_k_W, kmom);
        makeWallSink(application, q_pil_p, q_pil_p_W, pmom);*/

        /***********************************************************************
         * CONTRACTIONS: pi and K 2pt contractions with mom = p, k.
         **********************************************************************/
        // Wall-Point
        std::string PW_K_k = INIT_INDEX("PW_K_k", tK);
        std::string PW_K_p = INIT_INDEX("PW_K_p", tK);
        std::string PW_pi_k = INIT_INDEX("PW_pi_k", tpi);
        std::string PW_pi_p = INIT_INDEX("PW_pi_p", tpi);
        mesonContraction(application, 2, q_Kl_0, q_Ks_k, PW_K_k, kmom);
        mesonContraction(application, 2, q_Kl_0, q_Ks_p, PW_K_p, pmom);
        mesonContraction(application, 2, q_pil_k, q_pil_0, PW_pi_k, kmom);
        mesonContraction(application, 2, q_pil_p, q_pil_0, PW_pi_p, pmom);
        // Wall-Wall, to be done - requires modification of meson module.

        /***********************************************************************
         * CONTRACTIONS: 3pt Weak Hamiltonian, C & W (non-Eye type) classes.
         **********************************************************************/
        std::string HW_CW_k = LABEL_3PT("HW_CW_k", tK, tpi);
        std::string HW_CW_p = LABEL_3PT("HW_CW_p", tK, tpi);
        weakContractionNonEye(application, 3, q_Kl_0, q_Ks_k, q_pil_k, q_pil_0, HW_CW_k);
        weakContractionNonEye(application, 3, q_Kl_0, q_Ks_p, q_pil_p, q_pil_0, HW_CW_p);

        /***********************************************************************
         * CONTRACTIONS: 3pt sd insertion.
         **********************************************************************/

        for (unsigned int nn = 0; nn < n_noise; ++nn)
        {
            /*******************************************************************
             * CONTRACTIONS: 3pt Weak Hamiltonian, S and E (Eye type) classes.
             ******************************************************************/

        }

        // Perform separate contractions for each t_J position.
        for (unsigned int i = 0; i < tJs.size(); ++i)
        {
            // Sequential sources for current insertions. Local for now,
            // gamma_0 only.
            unsigned int tJ = (tJs[i] + tK) % nt;
            MSource::SeqGamma::Par seqPar;
            std::string q_KlCl_q   = LABEL_3PT("Q_KlCl_q", tK, tJ);
            std::string q_KsCs_mq  = LABEL_3PT("Q_KsCs_mq", tK, tJ);
            std::string q_pilCl_q  = LABEL_3PT("Q_pilCl_q", tpi, tJ);
            std::string q_pilCl_mq = LABEL_3PT("Q_pilCl_mq", tpi, tJ);
            MAKE_SEQUENTIAL_PROP(tJ, q_Kl_0, qmom, q_KlCl_q, solvers[light]);
            MAKE_SEQUENTIAL_PROP(tJ, q_Ks_k, mqmom, q_KsCs_mq, solvers[strange]);
            MAKE_SEQUENTIAL_PROP(tJ, q_pil_p, qmom, q_pilCl_q, solvers[light]);
            MAKE_SEQUENTIAL_PROP(tJ, q_pil_0, mqmom, q_pilCl_mq, solvers[light]);

            /*******************************************************************
             * CONTRACTIONS: pi and K 3pt contractions with current insertion.
             ******************************************************************/
            // Wall-Point
            std::string C_PW_Kl   = LABEL_3PT("C_PW_Kl", tK, tJ);
            std::string C_PW_Ksb  = LABEL_3PT("C_PW_Ksb", tK, tJ);
            std::string C_PW_pilb = LABEL_3PT("C_PW_pilb", tK, tJ);
            std::string C_PW_pil  = LABEL_3PT("C_PW_pil", tK, tJ);
            mesonContraction(application, 3, q_KlCl_q, q_Ks_k, C_PW_Kl, pmom);
            mesonContraction(application, 3, q_Kl_0, q_KsCs_mq, C_PW_Ksb, pmom);
            mesonContraction(application, 3, q_pil_0, q_pilCl_q, C_PW_pilb, kmom);
            mesonContraction(application, 3, q_pilCl_mq, q_pil_p, C_PW_pil, kmom);
            // Wall-Wall, to be done.

            /*******************************************************************
             * CONTRACTIONS: 4pt contractions, C & W classes.
             ******************************************************************/
            std::string CW_Kl   = LABEL_4PT("CW_Kl", tK, tJ, tpi);
            std::string CW_Ksb  = LABEL_4PT("CW_Ksb", tK, tJ, tpi);
            std::string CW_pilb = LABEL_4PT("CW_pilb", tK, tJ, tpi);
            std::string CW_pil  = LABEL_4PT("CW_pil", tK, tJ, tpi);
            weakContractionNonEye(application, 4, q_KlCl_q, q_Ks_k, q_pil_p, q_pil_0, CW_Kl);
            weakContractionNonEye(application, 4, q_Kl_0, q_KsCs_mq, q_pil_p, q_pil_0, CW_Ksb);
            weakContractionNonEye(application, 4, q_Kl_0, q_Ks_k, q_pilCl_q, q_pil_0, CW_pilb);
            weakContractionNonEye(application, 4, q_Kl_0, q_Ks_k, q_pil_p, q_pilCl_mq, CW_pil);

            /*******************************************************************
             * CONTRACTIONS: 4pt contractions, sd insertions.
             ******************************************************************/

            // Sequential sources for each noise propagator.
            for (unsigned int nn = 0; nn < n_noise; ++nn)
            {
                /***************************************************************
                 * CONTRACTIONS: 4pt contractions, S & E classes.
                 **************************************************************/

                /***************************************************************
                 * CONTRACTIONS: 4pt contractions, pi0 disconnected loop.
                 **************************************************************/
            }
        }
    }
    // execution
    application.saveParameterFile("rarekaon_000_100_tK0_tpi16_tJ8_noloop_mc0.2.xml");
    application.run();

    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
