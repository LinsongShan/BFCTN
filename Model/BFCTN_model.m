classdef BFCTN_model
    properties
        msi; hsi; Org; 
        rank; 
        P1; P2; P3;
        lambda_a_0; lambda_b_0;
        alpha_c_0; alpha_d_0;
        beta_e_0; beta_f_0; 
        W; w; H; h; S; s; N_k;
        dim; order; 
        Factors; lambda_12; lambda_13; lambda_14; lambda_23; lambda_24; lambda_34; alpha; beta;
        Lambda1; Lambda2; Lambda3; Lambda4;
        h1_1; h1_2; h2_1; h2_2; h3_1; h3_2; h4_1; h4_2; h5_1; h5_2; h6_1; h6_2;
        MHat; HHat; Hat; hat; init_hat;initRMSE; initPSNR; initerror;RMSE_List; PSNR_List;relative_error_map
        pre; cur; error_history
    end
    
    methods
        % create object
        function obj = BFCTN_model(msi, hsi, Org, P1, P2, P3, rank, hyperparameters)
            obj.msi = msi;
            obj.hsi = hsi;
            obj.Org = Org;
            obj.P1 = P1;
            obj.P2 = P2;
            obj.P3 = P3;
            obj.rank = rank;
            obj.lambda_a_0 = hyperparameters.lambda_a_0;
            obj.lambda_b_0 = hyperparameters.lambda_b_0;
            obj.alpha_c_0 = hyperparameters.alpha_c_0;
            obj.alpha_d_0 = hyperparameters.alpha_d_0;
            obj.beta_e_0 = hyperparameters.beta_e_0;
            obj.beta_f_0 = hyperparameters.beta_f_0;
        end

        % initialize parameter
        function self = initialize(self)
            self.W = size(self.msi, 1);
            self.w = size(self.hsi, 1);
            self.H = size(self.msi, 2);
            self.h = size(self.hsi, 2);
            self.S = size(self.hsi, 3);
            self.s = size(self.msi, 3);
            self.N_k = size(self.msi, 4);
            self.dim = [self.W, self.H, self.S, self.N_k];  
            self.order = 4;    

            % request facotr matrix memory
            self.Factors = cell(self.order, 1);   % cell save Factor Matrix
            self.Factors{1} = zeros(self.rank(1,2), self.rank(1,3), self.rank(1,4), self.dim(1));
            self.Factors{2} = zeros(self.rank(1,2), self.rank(2,3), self.rank(2,4), self.dim(2));
            self.Factors{3} = zeros(self.rank(1,3), self.rank(2,3), self.rank(3,4), self.dim(3));
            self.Factors{4} = zeros(self.rank(1,4), self.rank(2,4), self.rank(3,4), self.dim(4));

            % Initialise lambda
            self.lambda_12 = ones(self.rank(1,2), 1);
            self.lambda_13 = ones(self.rank(1,3), 1);
            self.lambda_14 = ones(self.rank(1,4), 1);
            self.lambda_23 = ones(self.rank(2,3), 1);
            self.lambda_24 = ones(self.rank(2,4), 1);
            self.lambda_34 = ones(self.rank(3,4), 1);

            % Facotrs Matrix 
            %T1
            T_1 = double(tenmat(self.Factors{1}, 4));
            tmp = kr(self.lambda_14,self.lambda_13, self.lambda_12);
            for i = 1:self.dim(1)
                    T_1(i,:) = mvnrnd(zeros(self.rank(1,2) * self.rank(1,3) * self.rank(1,4), 1), diag(tmp));
            end
            self.Factors{1} = reshape(T_1',  [self.rank(1,2), self.rank(1,3), self.rank(1,4), self.dim(1)]);

            %T2
            T_2 = double(tenmat(self.Factors{2}, 4));
            for i = 1:self.dim(2)
                    T_2(i,:) = mvnrnd(zeros(self.rank(1,2) * self.rank(2,3) * self.rank(2,4), 1), diag(kr(self.lambda_24,self.lambda_23, self.lambda_12)));
            end
            self.Factors{2} = reshape(T_2',  [self.rank(1,2), self.rank(2,3), self.rank(2,4), self.dim(2)]);

            %T3
            T_3 = double(tenmat(self.Factors{3}, 4));
            for i = 1:self.dim(3)
                    T_3(i,:) = mvnrnd(zeros(self.rank(1,3) * self.rank(2,3) * self.rank(3,4), 1), diag(kr(self.lambda_34,self.lambda_23, self.lambda_13)));
            end
            self.Factors{3} = reshape(T_3',  [self.rank(1,3), self.rank(2,3), self.rank(3,4), self.dim(3)]);

            %T4
            T_4 = double(tenmat(self.Factors{4}, 4));
            for i = 1:self.dim(4)
                    T_4(i,:) = mvnrnd(zeros(self.rank(1,4) * self.rank(2,4) * self.rank(3,4), 1), diag(kr(self.lambda_34,self.lambda_24, self.lambda_14)));
            end
            self.Factors{4} = reshape(T_4',  [self.rank(1,4), self.rank(2,4), self.rank(3,4), self.dim(4)]);

            % Initialise alpha, beta
            self.alpha = 1;
            self.beta = 1;
        end

        %% run model
        function self = run(self, RUN_MAX_iterations)
            self.RMSE_List = zeros(RUN_MAX_iterations,1);
            self.PSNR_List = zeros(RUN_MAX_iterations,1);
            self.relative_error_map = zeros(RUN_MAX_iterations,1);
            self.init_hat = double(tnprod_new(self.Factors));
            self.initRMSE = sqrt( sum((self.Org(:)-self.init_hat(:)).^2)./length(self.Org(:)) );
            self.initPSNR = lyPSNR(self.Org, self.init_hat);
            self.initerror = calculateRelativeError3D(self.Org, self.init_hat);
            fprintf('                   init:   rmse:%g  psnr:%g \n',  self.initRMSE,self.initPSNR);
                      
            for iter=1:RUN_MAX_iterations
                % update Facotrs{1} - T_1
                self.Factors{1} = self.upgrade_T1();

                % update Facotrs{2} - T_2
                self.Factors{2} = self.upgrade_T2();

                % update Facotrs{3} - T_3
                self.Factors{3} = self.upgrade_T3();

                % update Facotrs{4} - T_4
                self.Factors{4} = self.upgrade_T4();
                
                % update lambda
                self.lambda_12 = self.upgrade_lambda_a_12() ./ self.upgrade_lambda_b_12();
                self.lambda_13 = self.upgrade_lambda_a_13() ./ self.upgrade_lambda_b_13();
                self.lambda_14 = self.upgrade_lambda_a_14() ./ self.upgrade_lambda_b_14();
                self.lambda_23 = self.upgrade_lambda_a_23() ./ self.upgrade_lambda_b_23();
                self.lambda_24 = self.upgrade_lambda_a_24() ./ self.upgrade_lambda_b_24();
                self.lambda_34 = self.upgrade_lambda_a_34() ./ self.upgrade_lambda_b_34();

                % update alpha
                self.alpha =  self.upgrade_alpha_c() / self.upgrade_alpha_d();

                % alpha1 = self.alpha;

                % update beta
                self.beta = self.upgrade_beta_e() / self.upgrade_beta_f();

                
                %% update Z^hat
                self.hat = double(tnprod_new(self.Factors));
                self.RMSE_List(iter,1) = sqrt( sum((self.Org(:)-self.hat(:)).^2)./length(self.Org(:)) );
                self.PSNR_List(iter,1) = lyPSNR(self.Org, self.hat);
                if mod(iter, 1)==0
                    fprintf('                   vb: Iter%3d.  rmse:%g  psnr:%g  \n', iter, self.RMSE_List(iter),self.PSNR_List(iter));
                end

            end
            fprintf('                    rmse:%g  psnr:%g\n ',  self.RMSE_List(iter),self.PSNR_List(iter));
        end

        





        %% update posterior parameter
        % lambda_12_a
        function lambda_a = upgrade_lambda_a_12(self)
            lambda_a = self.lambda_a_0 + 0.5*(self.W*self.rank(1,3)*self.rank(1,4) + self.H * self.rank(2,3) * self.rank(2,4)) ;
            lambda_a = lambda_a .* ones(self.rank(1,2), 1);
        end

        % lambda_12_b
        function lambda_b = upgrade_lambda_b_12(self)
            self.h1_1 = kr(self.lambda_14, self.lambda_13, ones(self.W,1));
            self.h1_2 = kr(self.lambda_24, self.lambda_23, ones(self.H,1));
            tmp1 = double(tenmat(self.Factors{1}, 1)) .* double(tenmat(self.Factors{1}, 1));
            tmp2 = double(tenmat(self.Factors{2}, 1)) .* double(tenmat(self.Factors{2}, 1));
            tmp3 = self.h1_1' * tmp1';
            tmp4 = self.h1_2' * tmp2';
            tmp = tmp3' + tmp4';
            lambda_b = self.lambda_b_0 + 0.5*tmp;
        end

        % lambda_13
        function lambda_a = upgrade_lambda_a_13(self)
            lambda_a = self.lambda_a_0 + 0.5*(self.W*self.rank(1,2)*self.rank(1,4) + self.S * self.rank(2,3) * self.rank(3,4)) ;
            lambda_a = lambda_a .* ones(self.rank(1,3), 1);
        end

        function lambda_b = upgrade_lambda_b_13(self)
            self.h2_1 = kr(self.lambda_14, self.lambda_12, ones(self.W,1));
            self.h2_2 = kr(self.lambda_34, self.lambda_23, ones(self.S,1));

            tmp1 = double(tenmat(self.Factors{1}, 2)) .* double(tenmat(self.Factors{1}, 2));
            tmp2 = double(tenmat(self.Factors{3}, 1)) .* double(tenmat(self.Factors{3}, 1));
            tmp3 = self.h2_1' * tmp1';
            tmp4 = self.h2_2' * tmp2';
            tmp = tmp3' + tmp4';
            lambda_b = self.lambda_b_0 + 0.5*tmp;
        end

        % lambda_14
        function lambda_a = upgrade_lambda_a_14(self)
            lambda_a = self.lambda_a_0 + 0.5*(self.W*self.rank(1,2)*self.rank(1,3) + self.N_k * self.rank(2,4) * self.rank(3,4)) ;
            lambda_a = lambda_a .* ones(self.rank(1,4), 1);
        end

        function lambda_b = upgrade_lambda_b_14(self)
            self.h3_1 = kr(self.lambda_13, self.lambda_12, ones(self.W,1));
            self.h3_2 = kr(self.lambda_34, self.lambda_24, ones(self.N_k,1));

            tmp1 = double(tenmat(self.Factors{1}, 3)) .* double(tenmat(self.Factors{1}, 3));
            tmp2 = double(tenmat(self.Factors{4}, 1)) .* double(tenmat(self.Factors{4}, 1));
            tmp3 = self.h3_1' * tmp1';
            tmp4 = self.h3_2' * tmp2';
            tmp = tmp3' + tmp4';
            lambda_b = self.lambda_b_0 + 0.5*tmp;
        end

        % lambda_23
        function lambda_a = upgrade_lambda_a_23(self)
            lambda_a = self.lambda_a_0 + 0.5*(self.H*self.rank(1,2)*self.rank(2,4) + self.S * self.rank(1,3) * self.rank(3,4)) ;
            lambda_a = lambda_a .* ones(self.rank(2,3), 1);
        end

        function lambda_b = upgrade_lambda_b_23(self)
            self.h4_1 = kr(self.lambda_24, self.lambda_12, ones(self.H,1));
            self.h4_2 = kr(self.lambda_34, self.lambda_13, ones(self.S,1));

            tmp1 = double(tenmat(self.Factors{2}, 2)) .* double(tenmat(self.Factors{2}, 2));
            tmp2 = double(tenmat(self.Factors{3}, 2)) .* double(tenmat(self.Factors{3}, 2));
            tmp3 = self.h4_1' * tmp1';
            tmp4 = self.h4_2' * tmp2';
            tmp = tmp3' + tmp4';
            lambda_b = self.lambda_b_0 + 0.5*tmp;
        end

        % lambda_24
        function lambda_a = upgrade_lambda_a_24(self)
            lambda_a = self.lambda_a_0 + 0.5*( self.H * self.rank(1,2)*self.rank(2,3) + self.N_k * self.rank(1,4) * self.rank(3,4)) ;
            lambda_a = lambda_a .* ones(self.rank(2,4), 1);
        end

        function lambda_b = upgrade_lambda_b_24(self)
            self.h5_1 = kr(self.lambda_23, self.lambda_12, ones(self.H,1));
            self.h5_2 = kr(self.lambda_34, self.lambda_14, ones(self.N_k,1));

            tmp1 = double(tenmat(self.Factors{2}, 3)) .* double(tenmat(self.Factors{2}, 3));
            tmp2 = double(tenmat(self.Factors{4}, 2)) .* double(tenmat(self.Factors{4}, 2));
            tmp3 = self.h5_1' * tmp1';
            tmp4 = self.h5_2' * tmp2';
            tmp = tmp3' + tmp4';
            lambda_b = self.lambda_b_0 + 0.5*tmp;
        end

        % lambda_34
        function lambda_a = upgrade_lambda_a_34(self)
            lambda_a = self.lambda_a_0 + 0.5*(self.S * self.rank(1,3)*self.rank(2,3) + self.N_k * self.rank(1,4) * self.rank(2,4)) ;
            lambda_a = lambda_a .* ones(self.rank(3,4), 1);
        end

        function lambda_b = upgrade_lambda_b_34(self)
            self.h6_1 = kr(self.lambda_23, self.lambda_13, ones(self.S,1));
            self.h6_2 = kr(self.lambda_24, self.lambda_14, ones(self.N_k,1));

            tmp1 = double(tenmat(self.Factors{3}, 3)) .* double(tenmat(self.Factors{3}, 3));
            tmp2 = double(tenmat(self.Factors{4}, 3)) .* double(tenmat(self.Factors{4}, 3));
            tmp3 = self.h6_1' * tmp1';
            tmp4 = self.h6_2' * tmp2';
            tmp = tmp3' + tmp4';
            lambda_b = self.lambda_b_0 + 0.5*tmp;
        end

       
        % alpha_c
        function alpha_c = upgrade_alpha_c(self)
            alpha_c = self.alpha_c_0 + (self.W * self.H * self.s * self.N_k)/2;
        end
        
        % alpha_d
        function alpha_d = upgrade_alpha_d(self)
            M_Factors = self.Factors;
            M_Factors{3}  = ttm(tensor(M_Factors{3}), self.P3, 4);
            M_Factors{3} = double(M_Factors{3});

            temp = 0.5*((self.msi-double(tnprod_new(M_Factors))).^2);
            alpha_d = self.alpha_d_0 + sum(temp(:));
        end

        % beta_e
        function beta_e = upgrade_beta_e(self)
            beta_e = self.beta_e_0 + (self.w * self.h * self.S * self.N_k)/2;
        end
        
        % beta_f
        function beta_f = upgrade_beta_f(self)
            H_Factors = self.Factors;
            H_Factors{1}  = ttm(tensor(H_Factors{1}), self.P1, 4);
            H_Factors{1} = double(H_Factors{1});
            H_Factors{2}  = ttm(tensor(H_Factors{2}), self.P2, 4);
            H_Factors{2} = double(H_Factors{2});

            temp = 0.5*( (self.hsi-double(tnprod_new(H_Factors))).^2);
            beta_f = self.beta_f_0 + sum(temp(:));
        end


        % Factor T1
        function res = upgrade_T1(self)
            M_Factors = self.Factors;
            M_Factors{3}  = ttm(tensor(M_Factors{3}), self.P3, 4);
            M_Factors{3} = double(M_Factors{3});

            H_Factors = self.Factors;
            H_Factors{1}  = ttm(tensor(H_Factors{1}), self.P1, 4);
            H_Factors{1} = double(H_Factors{1});
            H_Factors{2}  = ttm(tensor(H_Factors{2}), self.P2, 4);
            H_Factors{2} = double(H_Factors{2});

            self.Lambda1 = kr(self.lambda_14, self.lambda_13, self.lambda_12);

            var1 = tnreshape_new(tnprod_rest_new(H_Factors,1),4);
            var2 = tnreshape_new(tnprod_rest_new(M_Factors,1),4);
            
            tmp1 =  self.beta * self.P1' * self.P1;
            tmp2 = var1* var1';
            tmp3 = self.alpha * (var2 * var2') + diag(self.Lambda1);
            tmp4 = self.beta * self.P1' * double(tenmat(self.hsi,1))* var1' + self.alpha * double(tenmat(self.msi, 1)) * var2';

            % Solve Sylvester equations  :   tmp1*X*tmp2  +  X*tmp3 = tmp4  --->   tmpa*X + X*tmpb = tmpc; 
            tmpa = tmp1;
            tmpb = tmp3 / tmp2;
            tmpc = tmp4 / tmp2;
            res = sylvester(tmpa, tmpb, tmpc);

            res = reshape(res', [self.rank(1,2), self.rank(1,3), self.rank(1,4), self.dim(1)]);
        end

        % Factor T2
        function res = upgrade_T2(self)
            M_Factors = self.Factors;
            M_Factors{3}  = ttm(tensor(M_Factors{3}), self.P3, 4);
            M_Factors{3} = double(M_Factors{3});

            H_Factors = self.Factors;
            H_Factors{1}  = ttm(tensor(H_Factors{1}), self.P1, 4);
            H_Factors{1} = double(H_Factors{1});
            H_Factors{2}  = ttm(tensor(H_Factors{2}), self.P2, 4);
            H_Factors{2} = double(H_Factors{2});

            self.Lambda2 = kr(self.lambda_24, self.lambda_23, self.lambda_12);

            var1 = tnreshape_new(tnprod_rest_new(H_Factors,2),4);
            var2 = tnreshape_new(tnprod_rest_new(M_Factors,2),4);
     
            tmp1 =  self.beta * self.P2' * self.P2;
            tmp2 = var1 * var1';
            tmp3 = self.alpha * (var2 * var2') + diag(self.Lambda2);
            tmp4 = self.beta * self.P1' * double(tenmat(self.hsi,2))* var1' + self.alpha * double(tenmat(self.msi, 2)) * var2';

            % Solve Sylvester equations  :   tmp1*X*tmp2  +  X*tmp3 = tmp4  --->   tmpa*X + X*tmpb = tmpc; 
            tmpa = tmp1;
            tmpb = tmp3 / tmp2;
            tmpc = tmp4 / tmp2;
            res = sylvester(tmpa, tmpb, tmpc);

            res = reshape(res', [self.rank(1,2), self.rank(2,3), self.rank(2,4), self.dim(2)]);
        end

        % Factor T3
        function res = upgrade_T3(self)
            M_Factors = self.Factors;
            M_Factors{3}  = ttm(tensor(M_Factors{3}), self.P3, 4);
            M_Factors{3} = double(M_Factors{3});

            H_Factors = self.Factors;
            H_Factors{1}  = ttm(tensor(H_Factors{1}), self.P1, 4);
            H_Factors{1} = double(H_Factors{1});
            H_Factors{2}  = ttm(tensor(H_Factors{2}), self.P2, 4);
            H_Factors{2} = double(H_Factors{2});

            self.Lambda3 = kr(self.lambda_34, self.lambda_23, self.lambda_13);

            var1 = tnreshape_new(tnprod_rest_new(M_Factors,3),4);
            var2 = tnreshape_new(tnprod_rest_new(H_Factors,3),4);
            
            tmp1 =  self.alpha * self.P3' * self.P3;
            tmp2 = var1 * var1';
            tmp3 = self.beta * (var2 * var2') + diag(self.Lambda3);
            tmp4 = self.alpha * self.P3' * double(tenmat(self.msi,3))* var1' + self.beta * double(tenmat(self.hsi, 3)) * var2';

            % Solve Sylvester equations  :   tmp1*X*tmp2  +  X*tmp3 = tmp4  --->   tmpa*X + X*tmpb = tmpc; 
            tmpa = tmp1;
            tmpb = tmp3 / tmp2;
            tmpc = tmp4 / tmp2;
            res = sylvester(tmpa, tmpb, tmpc);

            res = reshape(res', [self.rank(1,3), self.rank(2,3), self.rank(3,4), self.dim(3)]);
        end


         % Factor T4
        function res = upgrade_T4(self)
            M_Factors = self.Factors;
            M_Factors{3}  = ttm(tensor(M_Factors{3}), self.P3, 4);
            M_Factors{3} = double(M_Factors{3});

            H_Factors = self.Factors;
            H_Factors{1}  = ttm(tensor(H_Factors{1}), self.P1, 4);
            H_Factors{1} = double(H_Factors{1});
            H_Factors{2}  = ttm(tensor(H_Factors{2}), self.P2, 4);
            H_Factors{2} = double(H_Factors{2});

            self.Lambda4 = kr(self.lambda_34, self.lambda_24, self.lambda_14);

            var1 = tnreshape_new(tnprod_rest_new(M_Factors,4),4);
            var2 = tnreshape_new(tnprod_rest_new(H_Factors,4),4);
            
            tmp1 =  self.alpha * (var1 * var1');
            tmp2 = self.beta * (var2 * var2') + diag(self.Lambda4);
            tmp3 = self.alpha * double(tenmat(self.msi, 4)) * var1' + self.beta * double(tenmat(self.hsi, 4)) * var2';
            
            res = tmp3 / (tmp1 + tmp2);

            res = reshape(res',[self.rank(1,4), self.rank(2,4), self.rank(3,4), self.dim(4)]);
        end


    end
end


