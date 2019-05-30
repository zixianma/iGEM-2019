
d3.csv("https://raw.githubusercontent.com/zixianma/iGEM-2019/master/core_dataset.csv", function (data) {
    var new_data = filter(data);
    console.log(new_data.length);
    // var csv = generateSequenceCSV(new_data)

    // var blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    // var a = document.createElement('a');
    // a.download = 'AcrVIB.csv';
    // a.href = URL.createObjectURL(blob);
    // a.click();

    // var len_dict = getSeqLengthsDistAndMax(new_data)[0];
    // var max_len = getSeqLengthsDistAndMax(new_data)[1];
    // var target_len = generateLength(len_dict);
    // var generated_seq = generateAAAsIfDependent(new_data, target_len);
    // while (new_data.includes(generated_seq)){
    //     generated_seq = generateAAAsIfDependent(new_data, target_len);
    // };
});

function generateSequenceCSV(seq_list) {
    var new_seq_list = []
    seq_list.forEach(function (seq) {
        var new_seq = ""
        for (var i = 0; i < seq.length; i++) {
            new_seq = new_seq + seq[i] + "-";
        }
        new_seq += "end,1";
        // console.log(new_seq);
        new_seq_list.push(new_seq);
    });
    let csvContent = "";

    new_seq_list.forEach(function (row) {
        // let row = rowArray.join(",");
        csvContent += row + "\r\n";
    });
    return csvContent
}

function filter(data, family) {
    data_array = [];
    for (var i = 0; i < data.length; i++) {
        // Can add contional statements to filter the existing protein database
        acr = data[i];
        if (acr.Family == family) {
            var aa_seq = acr.Seq;
            data_array.push(aa_seq);
        }
    }
    return data_array;
}

function getSeqLengthsDistAndMax(seq_list) {
    var len_dict = {};
    var max_len = 0;
    var total = seq_list.length;
    seq_list.forEach(function (seq) {
        //Set the maximum length
        if (seq.length > max_len) {
            max_len = seq.length;
        }
        //Add to the probabiliy of a length 
        var key = seq.length;
        len_dict[key] = len_dict[key] == undefined ? 1 : (len_dict[key] + 1);
    });

    var len_dict_keys = Object.keys(len_dict)
    len_dict_keys.forEach(function (len) {
        len_dict[len] = len_dict[len] / total
    });
    // console.log(len_dict)
    return [len_dict, max_len];
}

function generateLength(len_dict) {
    total_aa = 0;
    prob = Math.random();
    p_len = 0.0;
    //For each possible length, get its probability and return a length according 
    //to the randomly generated number
    var len_dict_keys = Object.keys(len_dict);
    for (var i = 0; i < len_dict_keys.length; i++) {
        len_key = len_dict_keys[i]
        p_len += len_dict[len_key];
        if (prob <= p_len) {
            total_aa = len_key;
            break;
        }
    }
    return total_aa;
}


function getNextDist(seq_list, index, index_aa) {
    //Filter the list of sequences by including only those with the specific a.a. at the 
    //specific index
    new_seq_list = [];
    if (index > 0) {
        seq_list.forEach(function (seq) {
            if (seq.length > index && seq[index - 1] == index_aa) {
                new_seq_list.push(seq);
            }
        });
    } else {
        new_seq_list = seq_list;
    };
    // console.log(new_seq_list)
    next_dist = {}

    new_seq_list.forEach(function (seq) {
        //Make sure that the starting index doesn't exceed the length of the sequence

        var key = seq[index];
        next_dist[key] = next_dist[key] == undefined ? 1 : (next_dist[key] + 1);
    });

    //Convert the counts distribution to a probability distribution
    var total = 0;
    var aa_keys = Object.keys(next_dist);
    aa_keys.forEach(function (aa) {
        total += next_dist[aa];
    });
    aa_keys.forEach(function (aa) {
        next_dist[aa] = next_dist[aa] / total;
    });

    return [new_seq_list, next_dist];
}


function generateAAAsIfDependent(seq_list, aa_length) {
    var result = "";
    console.log(aa_length)
    var continuous_prob_dist = {}
    for (var i = 0; i < aa_length; i++) {

        //Make sure that there's either something in "result" or there isn't a starting a.a.
        try {
            var prev_aa = i > 0 ? result[i - 1] : "";
        } catch (e) {
            if (e instanceof RangeError) {
                console.log("Oops! Dead end encountered");
            } else {
                console.log("Sorry, something is wrong.");
            }
        }
        //Get the distribution 
        var results = getNextDist(seq_list, i, prev_aa);
        seq_list = results[0]
        next_dist = results[1]
        console.log(next_dist)
        //Use the same scheme to add a.a. to the sequence
        var aa_keys = Object.keys(next_dist);
        var p = Math.random();
        var p_aa = 0.0;

        for (var j = 0; j < aa_keys.length; j++) {
            var aa_key = aa_keys[j]
            p_aa += next_dist[aa_key];
            if (prob <= p_aa) {
                result += aa_key;
                break;
            }
        }
    }
    continuous_prob_dist.name
    var key = result[i];
    next_dist[key] = next_dist[key] == undefined ? key : next_dist[key];
    continuous_prob_dist
    console.log(result)
    // return result;
}
