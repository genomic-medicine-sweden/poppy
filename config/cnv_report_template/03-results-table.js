class ResultsTable extends EventTarget {
  #table;
  #isFiltered;
  #header;
  #body;
  #data;
  #activeCaller;
  #tooltip;

  constructor(element, config) {
    super();

    if (!element) {
      throw Error("no d3 selection supplied");
    }

    this.#table = element;
    this.#isFiltered = config?.filter;

    this.#header = this.#table.append("thead").append("tr");
    this.#body = this.#table.append("tbody");

    if (!config?.data) {
      throw Error("no data supplied");
    }

    this.#data = config?.data;
    this.#activeCaller = config?.caller ? config.caller : 0;

    this.#tooltip = this.initTooltip();

    this.update();
  }

  set activeCaller(index) {
    this.#activeCaller = index;
    this.update();
  }

  set filter(isFiltered) {
    this.#isFiltered = isFiltered;
    this.update();
  }

  columnDef(col) {
    switch (col) {
      // Integers
      case "start":
      case "length":
      case "end":
        return {
          class: "right",
          format: (x) => x.toLocaleString(undefined, {}),
        };

      // Floating point, allow for missing numbers
      case "copyNumber":
      case "corrCopyNumber":
      case "baf":
      case "cn":
        return {
          class: "right tooltip-trigger",
          format: (x) => {
            if (x !== null && !isNaN(x)) {
              return x.toLocaleString(undefined, { minimumFractionDigits: 2 });
            }
            return "NA";
          },
        };

      case "view":
        return {
          class: "view-region-link",
          format: (x) => x,
        };

      case "type":
        return {
          class: "left",
          format: (x) => x,
        };

      case "genes":
        return {
          class: "left",
          format: (x) => x.join(", "),
        };

      // Strings
      case "caller":
      case "chromosome":
      default:
        return {
          class: "left",
          format: (x) => x,
        };
    }
  }

  columnLabel(col) {
    const columns = {
      view: "View",
      caller: "Caller",
      chromosome: "Chr",
      genes: "Genes",
      start: "Start",
      end: "End",
      length: "Length",
      type: "Type",
      cn: "CN",
      baf: "BAF",
    };

    if (columns[col]) {
      return columns[col];
    }

    return col;
  }

  initTooltip() {
    return d3
      .select("body")
      .append("div")
      .attr("class", "copy-number-tooltip hidden")
      .call((d) =>
        d
          .append("table")
          .call((t) =>
            t
              .append("thead")
              .append("tr")
              .selectAll("th")
              .data(["caller", "type", "cn"])
              .join("th")
              .text(this.columnLabel)
          )
          .call((t) => t.append("tbody"))
      );
  }

  copyNumberTooltip(data) {
    this.#tooltip.classed("hidden", false);
    this.#tooltip
      .select("tbody")
      .selectAll("tr")
      .data(data)
      .join("tr")
      .selectAll("td")
      .data((d) => Object.entries(d))
      .join("td")
      .text((d) => this.columnDef(d[0]).format(d[1]));
  }

  showCopyNumberTooltip(rowIndex, data) {
    this.copyNumberTooltip(data[rowIndex].others);
  }

  hideCopyNumberTooltip() {
    this.#tooltip.classed("hidden", true);
  }

  update() {
    let tableData = this.#data
      .map((d) =>
        d.callers[this.#activeCaller].cnvs
          .filter(
            (di) => !this.#isFiltered || (this.#isFiltered && di.passed_filter)
          )
          .map((di) => {
            const allCols = { view: "🔍", chromosome: d.chromosome, ...di };
            // Don't display caller and filter status in table
            const { caller, passed_filter, ...cols } = allCols;
            return cols;
          })
      )
      .flat();

    if (tableData.length === 0) {
      tableData = [{ "No data to display": [] }];
    } else {
      // Find the corresponding copy numbers from the other caller(s)
      for (let cnv of tableData) {
        // Same chromosome
        let chromData = this.#data.filter(
          (d) => d.chromosome === cnv.chromosome
        );
        // ... different caller
        let callerData = chromData[0].callers.filter(
          (_, i) => i !== this.#activeCaller
        );
        // ... same gene
        let otherCnvs = callerData
          .map((d) =>
            d.cnvs
              .filter(
                (c) => cnv.genes.filter((g) => c.genes.includes(g)).length > 0
              )
              .map((c) => {
                return { caller: d.name, type: c.type, cn: c.cn };
              })
          )
          .flat();
        cnv.others = otherCnvs;
      }
    }

    this.#header
      .selectAll("th")
      .data(
        Object.keys(tableData[0]).filter((k) => k !== "others"),
        (d) => d
      )
      .join("th")
      .text(this.columnLabel)
      .attr("class", (d) => this.columnDef(d).class);

    this.#body
      .selectAll("tr")
      .data(tableData)
      .join("tr")
      .attr("data-chromosome", (d) => d.chromosome)
      .attr("data-start", (d) => d.start)
      .attr("data-length", (d) => d.length)
      .attr("data-index", (_, i) => i)
      .selectAll("td")
      .data((d) => Object.entries(d).filter(([k, _]) => k !== "others"))
      .join("td")
      .text(([key, value]) => this.columnDef(key).format(value))
      .attr("class", ([key, _]) => this.columnDef(key).class);

    this.#body.selectAll(".view-region-link").on("click", (e) => {
      this.dispatchEvent(
        new CustomEvent("zoom-to-region", {
          detail: {
            chromosome: e.target.parentElement.dataset.chromosome,
            start: Number(e.target.parentElement.dataset.start),
            length: Number(e.target.parentElement.dataset.length),
          },
        })
      );
    });

    this.#body
      .selectAll(".tooltip-trigger")
      .on("mouseenter", (e) =>
        this.showCopyNumberTooltip(e.target.parentNode.dataset.index, tableData)
      )
      .on("mouseout", () => this.hideCopyNumberTooltip())
      .on("mousemove", (e) => {
        const dims = this.#tooltip.node().getBoundingClientRect();
        const offset = d3.select(".container").node().getBoundingClientRect().y;
        const maxHeight =
          document.documentElement.clientHeight - offset - dims.height;
        this.#tooltip
          .style("top", `${Math.min(e.layerY, maxHeight)}px`)
          .style("left", `${e.layerX - dims.width - 20}px`);
      });
  }
}
