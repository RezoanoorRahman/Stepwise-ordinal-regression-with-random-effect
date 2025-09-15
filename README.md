While working for a project, I wanted to use stepwise variable selection for an ordinal output in the presence of a random effect. 

*This is not a direct package;* it uses R package **ordinal** for fitting the model and **dplyr** for data manipulation.

For comparing multiple models, I used Akaike's information criterion (aic). 

Algorithm sketch:

- Start with the null model
- 
## Stepwise Selection Algorithm for Ordinal Regression with Fixed and Random Effects

**Input:**
- `Y`: ordinal response variable  
- `Fixed_Variables`: list of candidate fixed effect variables  
- `Random_Variables`: list of candidate random effect grouping variables  
- `Data`: dataset containing all variables  

---

**Algorithm:**

1. **Initialize**
   - `Old_Variables ← {Y}`
   - `Fixed_In ← ∅`
   - `Random_In ← ∅`
   - Fit baseline model:  
     `Old_Model ← clm(Y ~ 1, Data)`  
     `Old_AIC ← AIC(Old_Model)`

2. **Repeat until no improvement or no variables left**
   - **For each candidate in Fixed_Variables:**
     - `Select_Fixed ← Old_Variables ∪ {candidate}`
     - Fit model:  
       - If `Random_In = ∅`: `Fixed_Model ← clm(Y ~ Fixed_In ∪ candidate, Data)`  
       - Else: `Fixed_Model ← clmm(Y ~ Fixed_In ∪ candidate + (Random_In), Data)`
     - Compute improvement:  
       `Fixed_Improvement ← Old_AIC – AIC(Fixed_Model)`
   - Record best fixed candidate and its improvement.

   - **For each candidate in Random_Variables:**
     - `Select_Random ← Old_Variables ∪ {candidate}`
     - Fit model:  
       `Random_Model ← clmm(Y ~ Fixed_In + (Random_In ∪ candidate), Data)`
     - Compute improvement:  
       `Random_Improvement ← Old_AIC – AIC(Random_Model)`
   - Record best random candidate and its improvement.

   - **Compare best fixed vs. best random:**
     - `Best_Overall_Improvement ← max(Best_Fixed_Improvement, Best_Random_Improvement)`

   - **If Best_Overall_Improvement > 0:**
     - If `Best_Fixed_Improvement ≥ Best_Random_Improvement`:  
       - `Old_Model ← Best_Fixed_Model`  
       - `Fixed_In ← Fixed_In ∪ {Best_Fixed_Variable}`  
       - Remove `Best_Fixed_Variable` from `Fixed_Variables`
     - Else:  
       - `Old_Model ← Best_Random_Model`  
       - `Random_In ← Random_In ∪ {Best_Random_Variable}`  
       - Remove `Best_Random_Variable` from `Random_Variables`
     - Update:  
       `Old_AIC ← AIC(Old_Model)`  
       `Old_Variables ← Old_Variables ∪ {Selected_Variable}`

   - **Else:**  
     - Stop loop.

3. **Output**
   - `Final_Model ← Old_Model`
   - `Selected_Fixed ← Fixed_In`
   - `Selected_Random ← Random_In`
   - `Final_AIC ← Old_AIC`

---

**Notes:**
- `clm` is used when no random effects are present.  
- `clmm` is used once at least one random effect is included.  
- Greedy slection: at each step, choose the variable (fixed or random) with the largest reduction in AIC.  
