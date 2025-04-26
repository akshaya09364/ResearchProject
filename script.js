// script.js - Interactive Functionality for scRNA-seq Reanalysis Website

/**
 * Main Application Controller
 */
class ScRNAApp {
    constructor() {
        // DOM Elements
        this.dom = {
            geneSearch: document.getElementById('geneSearch'),
            searchBtn: document.getElementById('searchBtn'),
            geneResult: document.getElementById('geneResult'),
            umapViz: document.getElementById('umapViz'),
            zoomIn: document.getElementById('zoomIn'),
            zoomOut: document.getElementById('zoomOut'),
            resetView: document.getElementById('resetView')
        };

        // Data stores
        this.data = {
            genes: null,
            umap: null,
            currentZoom: 1,
            clusterInfo: {
                0: { name: "Neurons", color: '#4285f4', markerGenes: ["SNAP25", "SYT1"] },
                1: { name: "Astrocytes", color: '#34a853', markerGenes: ["GFAP", "AQP4"] },
                2: { name: "Oligodendrocytes", color: '#ea4335', markerGenes: ["OLIG2", "MOG"] },
                3: { name: "Microglia", color: '#fbbc05', markerGenes: ["CX3CR1", "P2RY12"] },
                4: { name: "Endothelial", color: '#673ab7', markerGenes: ["CLDN5", "FLT1"] },
                5: { name: "OPCs", color: '#ff5722', markerGenes: ["PDGFRA", "CSPG4"] }
            }
        };

        // Visualization constants
        this.constants = {
            pointRadius: 4,
            highlightRadius: 6,
            umapMargin: { top: 20, right: 20, bottom: 40, left: 40 },
            tooltip: null
        };

        // Initialize app
        this.init();
    }

    /**
     * Initialize application
     */
    async init() {
        // Load data
        try {
            await this.loadData();
            this.setupEventListeners();
            this.initTooltip();
            this.renderUMAP();
        } catch (error) {
            console.error('Initialization error:', error);
            this.showError('Failed to initialize application. Please try again later.');
        }
    }

    /**
     * Initialize tooltip
     */
    initTooltip() {
        this.constants.tooltip = d3.select("body").append("div")
            .attr("class", "umap-tooltip")
            .style("opacity", 0);
    }

    /**
     * Load required data
     */
    async loadData() {
        // Gene information (marker genes for different cell types)
        this.data.genes = {
            "SNAP25": {
                expression: "High in neuronal clusters (0)",
                meanExpression: 8.2,
                function: "Synaptic vesicle fusion",
                markerFor: "Neurons",
                clusters: [0]
            },
            "GFAP": {
                expression: "Astrocyte-specific (cluster 1)",
                meanExpression: 6.7,
                function: "Astrocyte cytoskeleton",
                markerFor: "Astrocytes",
                clusters: [1]
            },
            "OLIG2": {
                expression: "Oligodendrocyte-specific (cluster 2)",
                meanExpression: 5.9,
                function: "Oligodendrocyte differentiation",
                markerFor: "Oligodendrocytes",
                clusters: [2]
            },
            "MOG": {
                expression: "Oligodendrocyte marker (cluster 2)",
                meanExpression: 4.8,
                function: "Myelin formation",
                markerFor: "Oligodendrocytes",
                clusters: [2]
            },
            "NEUROD6": {
                expression: "Neuronal marker (cluster 0)",
                meanExpression: 7.1,
                function: "Neuronal differentiation",
                markerFor: "Neurons",
                clusters: [0]
            },
            "AQP4": {
                expression: "Astrocyte marker (cluster 1)",
                meanExpression: 6.3,
                function: "Water channel protein",
                markerFor: "Astrocytes",
                clusters: [1]
            }
        };

        try {
            // Load real UMAP data from JSON file
            const response = await fetch('data/umap_clusters.json');
            if (!response.ok) {
                throw new Error(`Failed to load UMAP data: ${response.status}`);
            }
            this.data.umap = await response.json();

            // Process the loaded data
            this.data.umap = this.data.umap.map(d => ({
                ...d,
                cluster: parseInt(d.cluster),
                cellId: d.cell_id || `cell_${Math.random().toString(36).substr(2, 8)}`
            }));

            console.log(`Loaded ${this.data.umap.length} cells with clusters:`, 
                [...new Set(this.data.umap.map(d => d.cluster))].sort());

        } catch (error) {
            console.error('Error loading UMAP data:', error);
            throw new Error('Failed to load UMAP data. See console for details.');
        }
    }

    /**
     * Set up event listeners
     */
    setupEventListeners() {
        // Gene search functionality
        this.dom.searchBtn.addEventListener('click', () => this.handleGeneSearch());
        this.dom.geneSearch.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') this.handleGeneSearch();
        });

        // UMAP controls
        this.dom.zoomIn.addEventListener('click', () => this.zoomUMAP(1.2));
        this.dom.zoomOut.addEventListener('click', () => this.zoomUMAP(0.8));
        this.dom.resetView.addEventListener('click', () => this.resetUMAPView());
    }

    /**
     * Handle gene search
     */
    handleGeneSearch() {
        const gene = this.dom.geneSearch.value.trim().toUpperCase();
        if (!gene) {
            this.showError('Please enter a gene symbol');
            return;
        }

        if (this.data.genes[gene]) {
            this.displayGeneInfo(gene);
        } else {
            this.showError(`Gene "${gene}" not found. Try: ${Object.keys(this.data.genes).join(', ')}`);
        }
    }

    /**
     * Display gene information
     */
    displayGeneInfo(gene) {
        const info = this.data.genes[gene];
        const clusterInfo = info.clusters.map(c => 
            `Cluster ${c} (${this.data.clusterInfo[c]?.name || 'Unknown'})`
        ).join(', ');

        this.dom.geneResult.innerHTML = `
            <div class="gene-header">
                <h4>${gene}</h4>
                <span class="badge">Marker for: ${info.markerFor}</span>
            </div>
            <div class="gene-details">
                <p><strong>Expression Pattern:</strong> ${info.expression}</p>
                <p><strong>Expressed in:</strong> ${clusterInfo}</p>
                <p><strong>Mean Expression:</strong> ${info.meanExpression.toFixed(2)}</p>
                <p><strong>Function:</strong> ${info.function}</p>
            </div>
        `;
        
        // Highlight corresponding clusters in UMAP
        this.highlightClusters(gene);
    }

    /**
     * Highlight clusters where this gene is expressed
     */
    highlightClusters(gene) {
        const clustersToHighlight = this.data.genes[gene]?.clusters || [];
        
        d3.selectAll("#umap-svg circle")
            .attr("opacity", d => clustersToHighlight.includes(d.cluster) ? 1 : 0.2)
            .attr("r", d => clustersToHighlight.includes(d.cluster) ? 
                this.constants.highlightRadius : this.constants.pointRadius);
    }

    /**
     * Render UMAP visualization
     */
    renderUMAP() {
        // Clear previous visualization
        this.dom.umapViz.innerHTML = '';
        
        // Set up dimensions
        const width = this.dom.umapViz.clientWidth;
        const height = 300;
        const margin = this.constants.umapMargin;
        
        // Create SVG
        const svg = d3.select("#umapViz")
            .append("svg")
            .attr("id", "umap-svg")
            .attr("width", width)
            .attr("height", height);
        
        // Create main group
        const g = svg.append("g")
            .attr("transform", `translate(${margin.left},${margin.top})`);
        
        // Calculate scales
        const xExtent = d3.extent(this.data.umap, d => d.x);
        const yExtent = d3.extent(this.data.umap, d => d.y);
        
        this.scales = {
            x: d3.scaleLinear()
                .domain([xExtent[0] * 1.1, xExtent[1] * 1.1])
                .range([0, width - margin.left - margin.right]),
            y: d3.scaleLinear()
                .domain([yExtent[0] * 1.1, yExtent[1] * 1.1])
                .range([height - margin.top - margin.bottom, 0])
        };
        
        // Draw points
        g.selectAll("circle")
            .data(this.data.umap)
            .enter()
            .append("circle")
            .attr("cx", d => this.scales.x(d.x))
            .attr("cy", d => this.scales.y(d.y))
            .attr("r", this.constants.pointRadius)
            .attr("fill", d => this.data.clusterInfo[d.cluster]?.color || '#cccccc')
            .attr("opacity", 0.7)
            .on("mouseover", (event, d) => this.handlePointHover(event, d))
            .on("mouseout", () => this.handlePointOut());
        
        // Add axes
        this.drawAxes(g, width, height, margin);
        
        // Add legend
        this.drawLegend(svg, width, margin);
    }

    /**
     * Draw UMAP axes
     */
    drawAxes(g, width, height, margin) {
        // X axis
        g.append("g")
            .attr("class", "axis axis-x")
            .attr("transform", `translate(0,${height - margin.top - margin.bottom})`)
            .call(d3.axisBottom(this.scales.x).tickSizeOuter(0))
            .append("text")
            .attr("class", "axis-label")
            .attr("x", (width - margin.left - margin.right) / 2)
            .attr("y", 30)
            .attr("fill", "#000")
            .text("UMAP 1");
        
        // Y axis
        g.append("g")
            .attr("class", "axis axis-y")
            .call(d3.axisLeft(this.scales.y))
            .append("text")
            .attr("class", "axis-label")
            .attr("transform", "rotate(-90)")
            .attr("y", -30)
            .attr("x", -(height - margin.top - margin.bottom) / 2)
            .attr("fill", "#000")
            .text("UMAP 2");
    }

    /**
     * Draw cluster legend
     */
    drawLegend(svg, width, margin) {
        const legend = svg.append("g")
            .attr("class", "umap-legend")
            .attr("transform", `translate(${width - margin.right - 120},${margin.top})`);
        
        // Get unique clusters from data
        const clusters = [...new Set(this.data.umap.map(d => d.cluster))].sort();
        
        // Add title
        legend.append("text")
            .attr("x", 0)
            .attr("y", -10)
            .text("Cell Clusters")
            .style("font-weight", "bold");
        
        // Add legend items
        clusters.forEach((cluster, i) => {
            const clusterName = this.data.clusterInfo[cluster]?.name || `Cluster ${cluster}`;
            const clusterColor = this.data.clusterInfo[cluster]?.color || '#cccccc';
            
            const legendItem = legend.append("g")
                .attr("transform", `translate(0,${i * 20})`);
            
            legendItem.append("circle")
                .attr("r", 6)
                .attr("fill", clusterColor);
            
            legendItem.append("text")
                .attr("x", 12)
                .attr("y", 5)
                .text(clusterName)
                .style("font-size", "12px");
        });
    }

    /**
     * Handle point hover
     */
    handlePointHover(event, d) {
        d3.select(event.target)
            .attr("r", this.constants.highlightRadius)
            .attr("opacity", 1);
        
        // Get cluster info
        const clusterName = this.data.clusterInfo[d.cluster]?.name || `Cluster ${d.cluster}`;
        const markerGenes = this.data.clusterInfo[d.cluster]?.markerGenes?.join(", ") || "Unknown";
        
        // Update tooltip
        this.constants.tooltip
            .html(`
                <div class="tooltip-header">${clusterName}</div>
                <div><strong>Cell ID:</strong> ${d.cellId}</div>
                <div><strong>Markers:</strong> ${markerGenes}</div>
            `)
            .style("left", `${event.pageX + 15}px`)
            .style("top", `${event.pageY + 15}px`)
            .style("opacity", 1);
    }

    /**
     * Handle mouse out from point
     */
    handlePointOut() {
        d3.selectAll("#umap-svg circle")
            .attr("r", this.constants.pointRadius)
            .attr("opacity", 0.7);
        
        this.constants.tooltip.style("opacity", 0);
    }

    /**
     * Zoom UMAP visualization
     */
    zoomUMAP(factor) {
        this.data.currentZoom *= factor;
        d3.select("#umap-svg g")
            .attr("transform", `translate(${this.constants.umapMargin.left},${this.constants.umapMargin.top}) scale(${this.data.currentZoom})`);
    }

    /**
     * Reset UMAP view
     */
    resetUMAPView() {
        this.data.currentZoom = 1;
        d3.select("#umap-svg g")
            .attr("transform", `translate(${this.constants.umapMargin.left},${this.constants.umapMargin.top}) scale(1)`);
    }

    /**
     * Show error message
     */
    showError(message) {
        this.dom.geneResult.innerHTML = `
            <div class="error-message">
                <i class="fas fa-exclamation-circle"></i>
                <span>${message}</span>
            </div>
        `;
    }
}

// Initialize application when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    new ScRNAApp();
});