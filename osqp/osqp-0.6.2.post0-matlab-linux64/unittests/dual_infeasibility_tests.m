classdef dual_infeasibility_tests < matlab.unittest.TestCase

    properties
        P
        q
        A
        u
        l
        m,
        n
        solver
        options
        tol
    end

    methods(TestMethodSetup)
        function setup_problem(testCase)

            % Set options
            testCase.options = struct;
            testCase.options.verbose = 0;
            testCase.options.rho = 0.01;
            testCase.options.eps_abs = 1e-04;
            testCase.options.eps_rel = 1e-04;
            testCase.options.eps_prim_inf = 1e-15;  % Focus only on dual infeasibility
            testCase.options.check_termination = 1;

            % Setup tolerance
            testCase.tol = 1e-04;

        end
    end

    methods (Test)
        function test_dual_infeasible_lp(testCase)


            % unbouned example
            testCase.P = [];
            testCase.q = [2; -1];
            testCase.A = speye(2);
            testCase.l = [0; 0];
            testCase.u = [];

            % Setup solver
            testCase.solver = osqp;
            testCase.solver.setup(testCase.P, testCase.q, ...
                testCase.A, testCase.l, testCase.u, testCase.options);

            % Solve with OSQP
            results = testCase.solver.solve();

            % Check if they are close
            testCase.verifyEqual(results.info.status_val, ...
                                 testCase.solver.constant('OSQP_DUAL_INFEASIBLE'), ...
                                 'AbsTol', testCase.tol)

        end

        function test_dual_infeasible_qp(testCase)


            % unbouned example
            testCase.P = sparse(diag([4; 0]));
            testCase.q = [0; 2];
            testCase.A = sparse([1 1; -1 1]);
            testCase.l = [];
            testCase.u = [2; 3];

            % Setup solver
            testCase.solver = osqp;
            testCase.solver.setup(testCase.P, testCase.q, ...
                testCase.A, testCase.l, testCase.u, testCase.options);

            % Solve with OSQP
            results = testCase.solver.solve();

            % Check if they are close
            testCase.verifyEqual(results.info.status_val, ...
                                 testCase.solver.constant('OSQP_DUAL_INFEASIBLE'), ...
                                 'AbsTol', testCase.tol)

        end

        function test_primal_dual_infeasible(testCase)
            testCase.P = sparse(2,2);
            testCase.q = [-1;-1];
            testCase.A = sparse([1 -1; -1 1;1 0;0 1]);
            testCase.l = [1; 1; 0; 0];
            testCase.u = Inf(4,1);

            % Setup solver
            testCase.solver = osqp;
            testCase.solver.setup(testCase.P, testCase.q, ...
                testCase.A, testCase.l, testCase.u, testCase.options);

            % Set warm starting points to avoid primal infeasibility detection at
            % first step
            x = 25*ones(2,1);
            y = -2*ones(4,1);
            testCase.solver.warm_start('x', x, 'y', y);

            % Solve with OSQP
            results = testCase.solver.solve();

            % Check that the solver status is correct
            testCase.verifyThat(results.info.status_val, ...
                                matlab.unittest.constraints.IsSubsetOf([ ...
                                    testCase.solver.constant('OSQP_PRIMAL_INFEASIBLE'), ...
                                    testCase.solver.constant('OSQP_DUAL_INFEASIBLE')]));

        end

    end

end
