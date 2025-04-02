import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error


# STYLING / PAGE LAYOUT
# Set up app page settings
st.set_page_config(page_title="Cell Survival Ratio Simulator - PHYS866", layout="wide")
if "run_sim" not in st.session_state:
    st.session_state["run_sim"] = False

# Set up the page title and description
st.title("PHYS866: Cell Survival Ratio Simulation")
st.markdown("This is a simulation for a University assignment and should not be used for any diagnostic or clinical purposes under any circumstances - for any questions please contact Sam.McCaughran@Gmail.com.")

# FIXED VALUES
# Tissue Alpha and Beta Values
tissue_parameters = {
    "Bladder":      {"alpha": 0.044, "beta": 0.003},
    "Breast":       {"alpha": 0.103, "beta": 0.026},
    "Cervix":       {"alpha": 0.130, "beta": 0.003},
    "CNS":          {"alpha": 0.081, "beta": 0.014},
    "Head and Neck":{"alpha": 0.059,  "beta": 0.003},
    "Liver":        {"alpha": 0.037,  "beta": 0.003},
    "Oesophagus":   {"alpha": 0.04,  "beta": 0.008},
    "Prostate":     {"alpha": 0.067, "beta": 0.024},
    "Rectum":       {"alpha": 0.313, "beta": 0.050},
    "Skin":         {"alpha": 0.005, "beta": 0.009},
    "Custom":       {"alpha": None, "beta": None}
}



# Variable for export data
export_data_csv = None

# Sidebar inputs for simulation parameters
st.sidebar.markdown(
    f"""
    <div style="text-align: left;">
        <img src="https://github.com/Sam-McCaughran/PHYS866_Part1/blob/main/test2.png?raw=true" style="width: 80%; max-width: 1080px;">
    </div>
    """,
    unsafe_allow_html=True
)
st.sidebar.header("Simulation Parameters")
simulate = st.sidebar.toggle("Run Simulation", value=True)
show_scatter = st.sidebar.toggle("Show Scatter Points", value=True)
selected_tissue = st.sidebar.selectbox(
    "Tissue Type",
    options=list(tissue_parameters.keys())
)
if selected_tissue == "Custom":
    alpha_value = st.sidebar.number_input("Alpha", value=0.1, step=0.001, format="%.4f")
    beta_value = st.sidebar.number_input("Beta", value=0.01, step=0.001, format="%.4f")
else:
    alpha_value = tissue_parameters[selected_tissue]["alpha"]
    beta_value = tissue_parameters[selected_tissue]["beta"]
    st.sidebar.markdown(f"**α = {alpha_value:.4f}**  \n**β = {beta_value:.4f}**")
num_experiments = st.sidebar.slider("Number of Experiments", min_value=0, max_value=15, value=0)
dose_range = st.sidebar.slider("Dose Range", 1, 50, value=25, step=1)
doses = np.linspace(0, dose_range, 100)
y_log_scale = st.sidebar.checkbox("Logarithmic Y-axis", value=True)

# SIDE BAR NOISE
st.sidebar.header("Noise")
add_counting_noise = st.sidebar.checkbox("Add Counting Noise", value=False)
initial_cells = st.sidebar.slider("Initial Cell Count", 100, 10000, value=1000, step=100)
add_dose_error = st.sidebar.checkbox("Add Dose Delivery Error", value=False)
dose_error_std = st.sidebar.slider("Dose Error Std Dev (%)", 0.0, 10.0, 0.0, step=0.1)
add_ab_variation = st.sidebar.checkbox("Add α/β Variation", value=False)
alpha_sd = st.sidebar.slider("Alpha Std Dev", 0.0000, 0.1, 0.000, step=0.0001)
beta_sd = st.sidebar.slider("Beta Std Dev", 0.0000, 0.1, 0.000, step=0.0001)
hypoxia_level = st.sidebar.selectbox("Tumour Oxygenation", ["Normoxic (OER: 1)", "Moderate Hypoxia (OER: 2)", "Severe Hypoxia (OER: 3)"])

# Map dropdown to values
oer_map = {
    "Normoxic (OER: 1)": 1.0,
    "Moderate Hypoxia (OER: 2)": 2.0,
    "Severe Hypoxia (OER: 3)": 3.0
}


# LOGIC
# Define the cell survival model using the LQ equation
def cell_radiation_model(dose, alpha, beta):
    """
    Compute cell survival fraction using the LQ model.
    
    Parameters:
        dose : scalar or np.array
            Radiation dose(s).
        alpha : float
            Parameter representing single-strand breaks.
        beta : float
            Parameter representing double-strand breaks.
    
    Returns:
        np.array: Survival fraction(s) computed from the LQ model.
    """
    return np.exp(-alpha * dose - beta * (dose ** 2))



# SIMULATION LOOP
if simulate:

    alpha = alpha_value
    beta = beta_value
    export_data = []
    oer_value = oer_map[hypoxia_level]
    
    # Compute the theoretical survival fraction
    theoretical_sf = cell_radiation_model(doses, alpha_value, beta_value)
    
    # Initialize a Plotly figure
    fig = go.Figure()
    # Run experiments by adding noise to the theoretical model
    for i in range(num_experiments):

        # Sample alpha and beta if variation is enabled
        if add_ab_variation:
            alpha = np.random.normal(alpha_value, alpha_sd)
            beta = np.random.normal(beta_value, beta_sd)
            alpha = max(alpha, 0)
            beta = max(beta, 0)
        else:
            alpha = alpha_value
            beta = beta_value

        # Apply hypoxia effect to alpha and beta values
        alpha /= oer_value
        beta /= oer_value

        # Compute the new survival fraction
        new_sf = cell_radiation_model(doses, alpha, beta)

        # Simulate cell colonies
        expected_counts = new_sf * initial_cells

        # Simulate cell counting noise
        if add_counting_noise:
            noisy_counts = np.random.normal(loc=expected_counts, scale=np.sqrt(expected_counts))
            # Avoid negative numbers
            noisy_counts = np.clip(noisy_counts, 0, initial_cells)
        else:
            noisy_counts = expected_counts

        # Get out experimental data
        new_experiment = noisy_counts / initial_cells

        # Apply dose delivery error
        if add_dose_error:
            current_doses = doses + np.random.normal(0, dose_error_std / 100 * doses)
        else:
            current_doses = doses

        # Variable for output fit parameters 
        fitted_params = [] 

            # Fit the LQ model to the noisy data
        try:
            popt, _ = curve_fit(cell_radiation_model, current_doses, new_experiment, p0=(alpha, beta), bounds=(0, np.inf))
            fitted_alpha, fitted_beta = popt
            fitted_params.append((fitted_alpha, fitted_beta))
            
            # Plot the best fit curve
            fitted_sf = cell_radiation_model(current_doses, fitted_alpha, fitted_beta)

            # Get out base r2 and mean square error
            r2 = r2_score(new_experiment, fitted_sf)
            mse = mean_squared_error(new_experiment, fitted_sf)

            # Get out r2 and mean square error relative to the theoretical curve
            r2_relative = r2_score(new_experiment, theoretical_sf)
            mse_relative = mean_squared_error(new_experiment, theoretical_sf)

            # Plot fit lines
            fig.add_trace(go.Scatter(
                x=current_doses,
                y=fitted_sf,
                mode="lines",
                name=f"Fit {i+1} (α={fitted_alpha:.4f}, β={fitted_beta:.4f} r²: {r2:.4f})"
            ))
            
        except RuntimeError:
            st.warning(f"Fit failed for experiment {i+1}")
            continue
        
        # Plot scatter marks
        if show_scatter:
            fig.add_trace(go.Scatter(
                x=current_doses,
                y=new_experiment,
                mode="markers",
                name=f"Experiment {i+1}"
            ))

        export_data.append({
        "Experiment": i + 1,
        "Tissue Type": selected_tissue,
        "Sampled Alpha": alpha,
        "Sampled Beta": beta,
        "Fitted Alpha": fitted_alpha,
        "Fitted Beta": fitted_beta,
        "Counting Noise": add_counting_noise,
        "Cell Count": initial_cells,
        "Dose Error (%)": dose_error_std if add_dose_error else 0.0,
        "A/B Variation": add_ab_variation,
        "Alpha Std Dev": alpha_sd,
        "B Std Dev": beta_sd,
        "OER": oer_value,
        "R² - Relative to Fit": r2,
        "MSE - Relative to Fit": mse,
        "R² - Relative to Theoretical Curve": r2_relative,
        "MSE - Relative Theoretical Curve": mse_relative
    })
    
    # Customize the layout of the figure
    fig.update_layout(
        title="Cell Survival Fraction Plot: Theoretical vs Experimental",
        xaxis_title="Radiation Dose (Gy)",
        yaxis_title="Survival Fraction",
        template="plotly_white",
        height=650,
        #width=2000,
        yaxis_type="log" if y_log_scale else "linear"
    )

    # Add the theoretical LQ model curve
    fig.add_trace(go.Scatter(
        x=doses, 
        y=theoretical_sf, 
        mode="lines", 
        name=f"Theoretical LQ Model (α={alpha_value:.4f}, β={beta_value:.4f})",
        line=dict(color="black", dash="dash")
    ))

    # Display Plotly charT
    st.plotly_chart(fig, use_container_width=True)
    st.markdown("---")
    if export_data:
        df_export = pd.DataFrame(export_data)
        export_data_csv = df_export.to_csv(index=False).encode("utf-8")
        # Results display
        st.subheader("Experiment Summary")
        st.dataframe(df_export[["Experiment", "Fitted Alpha", "Fitted Beta", "R² - Relative to Theoretical Curve", "MSE - Relative Theoretical Curve"]])

# EXPORT BUTTON
st.markdown("---")
if export_data_csv:
    st.download_button(
            label="Download Experiment Data as CSV",
            data=export_data_csv,
            file_name="cell_survival_experiments.csv",
            mime="text/csv",
        )
    

