pls <- function(X,y,p_value_threshold,correlation_threshold) {
                        # x dataframe
                        df = data.frame(X)
                        df['Y']= y
                        # correlation matrix
                        corr = cor(df)


                        #scale data
                        df_scaled = scale(df)
                        df_scaled = data.frame(df_scaled)

                        n_features = ncol(X)
                        corr_y = corr[n_features + 1,1:n_features]
                        corr_yabs = abs(corr_y)

                        # Empty list to store the selected features
                        selected_features = list()

                        for (i in 1:n_features){
                            if (corr_yabs[i] > correlation_threshold){
                                selected_features[[length(selected_features) + 1]] = i
                            }
                        }

                        sum1 = 0
                        sum2 = 0
                        for (i in selected_features){
                            sum1 = sum1 + corr[n_features + 1 , 1:n_features][i] * df_scaled[i]
                            sum2 = sum2 + corr_yabs[i]**2
                        }

                        T1 = (1/sqrt(sum2)) * sum1
                        T = data.frame(T1 = T1)
                        Y_scaled = df_scaled[n_features + 1]

                        j = 2
                        while (TRUE){
                            model = list()
                            for (i in 1:n_features){
                                X_ = cbind(df_scaled[i],T)
                                model[[i]] = lm(formula = unlist(Y_scaled) ~ unlist(T) + unlist(df_scaled[i]) )                    
                            }

                            col2 = list()
                            for (i in 1:n_features){
                                p_value = summary(model[[i]])$coefficients[,4][3]
                                if (p_value < p_value_threshold){
                                    col2[[length(col2) + 1]] = i
                                }
                            }

                            if (length(col2) == 0){
                                break
                            }

                            residu = list()

                            for (i in col2){
                                print(i)
                                model = lm(formula = unlist(df_scaled[i]) ~ unlist(T) )
                                residu[[length(residu) + 1]] = resid(model)
                            }

                            X1n = list()
                            for (r in residu){
                                X1n[[length(X1n) + 1]] = r/var(r)
                            }

                            model2 = list()
                            for (i in 1:length(X1n)){
                                dt = T
                                dt[1] = X1n[[i]]
                                model2[[i]] = lm(formula = unlist(Y_scaled) ~ unlist(dt) )
                            }

                            ceofs = list()
                            for (i in 1:length(model2)){
                                ceofs[[length(ceofs) + 1]] = summary(model2[[i]])$coefficients[2]
                            }

                            sum1 = 0
                            sum2 = 0
                            for (i in 1:length(residu)){
                                sum1 = sum1 + ceofs[[i]] * residu[[i]]
                                sum2 = sum2 + ceofs[[i]]**2
                            }

                            T[paste(c("T", 2), collapse = " ")] = (1/sqrt(sum2)) * sum1
                            j = j + 1

                    }
                    T['Y'] = y
                    return(T)
}
