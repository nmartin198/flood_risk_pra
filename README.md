# Probabilistic Risk Assessment (PRA) for Sustainable Water Resources Management: A Future Flood Inundation Example

This repository provides the source code and other files for the manuscript **Probabilistic Risk Assessment (PRA) for Sustainable Water Resources Management: A Future Flood Inundation Example**. The event tree which provides the scenario for risk assessment is shown below.
<br/>
<figure>
    <img src="/assets/Fig-04-Event_Tree.png"
         width="900"
         alt="Event Tree">
</figure>
<br/>

In this study, a PRA for future flood inundation is presented to examine solutions for the scenario where a regulatory compliant land use configuration is unsustainable when planning uncertainty considerations, related to future climate and weather, are applied. PRA is evaluation of the event tree shown above. 

- Future climate and weather are obtained from a weather attribution constrained projection of [Frio Basin weather parameter time series](https://github.com/nmartin198/wattrib_wg_frio).
    - [Initiating event summary](https://github.com/nmartin198/flood_risk_pra/tree/main/Frio_Synthetic_Weather_WG)
- For the stochasic obstruction height branch of the event tree, a Generalized Extreme Value distribution provides the variate.
    - [Stochastic obstruction height Jupyter notebook](https://github.com/nmartin198/flood_risk_pra/tree/main/Jupyter_Lab/Blockage_CDF_Formulation.ipynb) 
- Flood inundation is simulated with the [MOD_FreeSurf2D model](https://github.com/nmartin198/MOD_FreeSurf2D)


## Results

All results are available at [**PRA_EventTree_Analysis**](https://github.com/nmartin198/flood_risk_pra/tree/main/PRA_EventTree_Analysis).


## Author

* **Nick Martin** nick.martin@alumni.stanford.edu
<br/>

## License

This project is licensed under the GNU Affero General Public License v.3.0 - see the [LICENSE](LICENSE) file for details.
<br/>
