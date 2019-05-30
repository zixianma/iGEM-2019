
d3.csv("https://raw.githubusercontent.com/zixianma/iGEM-2019/master/core_dataset.csv").then(function(data) {
    console.log(data[0]);
  });

function main(data){
    for(var i = 0; i < data.length; i++){
    // Can add contional statements to filter the existing protein database
    acr = data[i]
    acr_type = acr[2]
    if (acr_type == "AcrF7"){
        aa_seq = acr[6]
}
    if (aa_seq != "Seq") {
        aa_seq_list.append(aa_seq)
    }


    //Get a len_dict storing the distribution of the a.a. lengths and get the max_length 
    len_dict, max_len = getSeqLengthsDistAndMax(aa_seq_list)
    console.log("There are a total of %d different amino acid sequences of this protein." % (len(aa_seq_list)))
    //print("The lengths distribution of the proteins is %s." % (len_dict))
    console.log("The maximum length of these aa sequences is %d." % (max_len))

    //Get the dictionary aa_dist which stores the distribution of a.a.s at each position
    aa_dist = getAADistrubution(aa_seq_list)

    //Convert the counts of a.a.s into probabilities
    aa_dist = convertCountIntoProb(aa_dist)
    new_data = rearrangeDataForHeatmap(aa_dist,"acrVIB")
    //Find out at which positions there are more than one possible a.a.s, and those are 
    //the positions of interest
    pos_of_interest =[]
    for (pos in aa_dist){
        if (len(aa_dist[pos]) > 1){
            pos_of_interest.append(pos)
        }
    }
    // print(len(pos_of_interest),pos_of_interest)
    //print(aa_dist)

    //Randomly generates a length based on the a.a. lenghts distribution 
    target_len = generateLength(len_dict)
    console.log("The randomly generated target length is %d." % (target_len))
    //Generate a novel a.a. sequence based on the distributions of different a.a.s
    //at different locations, assuming that they are INDEPENDENT of each other
    generated_aa_seq,length = generateAA(aa_dist,pos_of_interest,target_len)
    while (generated_aa_seq in aa_seq_list) {
        generated_aa_seq,length = generateAA(aa_dist,pos_of_interest,target_len)
    }
    console.log("Assume that the amino acids at different positions are independently distributed, we have: %s, length: %d" % (generated_aa_seq,length))

    //Generate a novel a.a. sequence based on the distributions of different a.a.s
    // at different locations, assuming that they are DEPENDENT of each other
    generated_aa_seq_de = generateAAAsIfDependent(aa_seq_list,target_len)
    while (generated_aa_seq in aa_seq_list) {
        generated_aa_seq_de  = generateAAAsIfDependent(aa_seq_list,target_len)
    }
    console.log("If the amino acids at different positions are dependent, we have: %s, length: %d" % (generated_aa_seq_de,len(generated_aa_seq)))
    }

    //This fucntion converts the counts distribution of a.a.s at different locs
    //into probability distribution 
function convertCountIntoProb(aa_dict){
        for (pos in aa_dict){
            total = 0
            for (aa in aa_dict[pos]){
                total += aa_dict[pos][aa]
            }
            for (aa in aa_dict[pos]){
                aa_dict[pos][aa] = aa_dict[pos][aa] / total 
            }
        }
        return aa_dict
}

function generateLength(len_dict){
    total_aa = 0
    prob = random.random()
    p_len = 0.0
    //For each possible length, get its probability and return a length according 
    //to the randomly generated number
    for (length in list(len_dict.keys())){
        p_len += len_dict[length]
        if (prob <= p_len){
            total_aa = length
            break
        }
    }
    return total_aa
}


//This is a helper function for generateAAAsIfDependent, and it's getting the 
//distribution of a.a.s at the next position given start index and a.a.
function getNextDist(seq_list,index,index_aa){
    //Filter the list of sequences by including only those with the specific a.a. at the 
    //specific index
        if (index > 0){
            seq_list = [seq for seq in seq_list if len(seq) > index and seq[index-1] == index_aa]
        }
        next_dist = {}
        for (seq in seq_list){
            //Make sure that the starting index doesn't exceed the length of the sequence
            if (seq[index] in next_dist){
                next_dist[seq[index]] += 1 
            } else {
                next_dist.setdefault(seq[index]) 
                next_dist[seq[index]] = 0
                next_dist[seq[index]] += 1 
            }
        }
        //Convert the counts distribution to a probability distribution
        total = 0
        for (aa in next_dist) {
            total += next_dist[aa]
        }
        for (aa in next_dist) {
            next_dist[aa] = next_dist[aa] / total 
        }
        //print(index,start_aa,next_dist)
        return seq_list,next_dist
}

//This function is for generating a.a. sequences assuming that the a.a.s are dependently
//distributed, which is most likely the case
function generateAAAsIfDependent(seq_list,aa_length){
    result = ""
    for (i in range(aa_length)){
        //Make sure that there's either something in "result" or there isn't a starting a.a.
        prev_aa = result[i-1] if i > 0 else ""
        //Get the distribution 
        seq_list,next_dist = getNextDist(seq_list,i,prev_aa)
        //Use the same scheme to add a.a. to the sequence
        aa_list = list(next_dist.keys())
        p = random.random()
        p_aa = 0.0
        for (aa in aa_list){
            p_aa += next_dist[aa]
            if (p <= p_aa){
                result += aa
                break
            }
        }
    return result
    }
}

function visualizeProcess(){

    // set the dimensions and margins of the graph
    var margin = {top: 80, right: 30, bottom: 30, left: 100},
      width = 1800 - margin.left - margin.right,
      height = 1000 - margin.top - margin.bottom;
    
    // append the svg object to the body of the page
    var svg = d3.select("#my_dataviz")
    .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")");
    
    //Read the data
    var e = document.getElementById('acr_id');
    var acr_id = e.options[e.selectedIndex].text;
    // var acr_id = document.getElementById('url_prefix').innerHTML;
    var url = "https://raw.githubusercontent.com/zixianma/iGEM-2019/master/heatmap_data/heatmap_data_" + acr_id + ".csv"
    
    d3.csv(url, function(data) {
    
      // Labels of row and columns -> unique identifier of the column called 'group' and 'variable'
      var myGroups = d3.map(data, function(d){return d.group;}).keys()
      var myVars = d3.map(data, function(d){return d.variable;}).keys()
    
      // Build X scales and axis:
      var x = d3.scaleBand()
        .range([ 0, width ])
        .domain(myGroups)
        .padding(0.05);
      svg.append("g")
        .style("font-size", 15)
        .attr("transform", "translate(0," + height + ")")
        .call(d3.axisBottom(x).tickSize(0))
        .select(".domain").remove()
    
      // Build Y scales and axis:
      var y = d3.scaleBand()
        .range([ height, 0 ])
        .domain(myVars)
        .padding(0.05);
      svg.append("g")
        .style("font-size", 15)
        .call(d3.axisLeft(y).tickSize(0))
        .select(".domain").remove()
    
      // Build color scale
      var myColor = d3.scaleSequential()
        .interpolator(d3.interpolateInferno)
        .domain([1,100])
    
      // create a tooltip
      var tooltip = d3.select("#my_dataviz")
        .append("div")
        .style("opacity", 0)
        .attr("class", "tooltip")
        .style("background-color", "white")
        .style("border", "solid")
        .style("border-width", "2px")
        .style("border-radius", "5px")
        .style("padding", "5px")
    
      // Three function that change the tooltip when user hover / move / leave a cell
      var mouseover = function(d) {
        tooltip
          .style("opacity", 1)
        d3.select(this)
          .style("stroke", "black")
          .style("opacity", 1)
      }
      var mousemove = function(d) {
        tooltip
          .html("The exact value of<br>this cell is: " + d.value)
          .style("left", (d3.mouse(this)[0] - 10) + "px")
          .style("top", (d3.mouse(this)[1] + 450) + "px")
      }
      var mouseleave = function(d) {
        tooltip
          .style("opacity", 0)
        d3.select(this)
          .style("stroke", "none")
          .style("opacity", 0.8)
      }
    
      // add the squares
      svg.selectAll()
        .data(data, function(d) {return d.group+':'+d.variable;})
        .enter()
        .append("rect")
          .attr("x", function(d) { return x(d.group) })
          .attr("y", function(d) { return y(d.variable) })
          .attr("rx", 4)
          .attr("ry", 4)
          .attr("width", x.bandwidth() )
          .attr("height", y.bandwidth() )
          .style("fill", function(d) { return myColor(90-d.value * 30)} )
          .style("stroke-width", 4)
          .style("stroke", "none")
          .style("opacity", 0.8)
        .on("mouseover", mouseover)
        .on("mousemove", mousemove)
        .on("mouseleave", mouseleave)
    })
    
    var title = "A heatmap of the amino acids composition of the " + acr_id + " proteins."
    // Add title to graph
    svg.append("text")
            .attr("x", 0)
            .attr("y", -50)
            .attr("text-anchor", "left")
            .style("font-size", "22px")
            .text(title);
    
    // Add subtitle to graph
    svg.append("text")
            .attr("x", 0)
            .attr("y", -20)
            .attr("text-anchor", "left")
            .style("font-size", "14px")
            .style("fill", "grey")
            .style("max-width", 400)
            .text("The x-axis represents the position of an amino acid sequence.\
            The y-axis is the one-letter code of amino acid.");
    } 