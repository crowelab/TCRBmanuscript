function query(v, d, j, cdr3len, count) {

var cursor=db.synthetics.aggregate(
[
{$match:
    {
        "v_family": v,
        "d_family": d,
        "j_family": j,
        "cdr3_aa_len": cdr3len
    }
},
{$sample: {size: count}},
{$group:
    {
        _id: {"v": "$v_family", "d":"$d_family", "j":"$j_family", "cdr3":"$cdr3_aa"},
        "count": {$sum:1}
    }
}
],
{allowDiskUse: true}
)

while(cursor.hasNext()) {
  var rec=cursor.next();
  print(rec._id.v + " " + rec._id.d + " " + rec._id.j + " " + rec._id.cdr3 + " " + rec.count);
}

}

var num_samples=1000;

for (var s=1; s<num_samples+1; s++) {
  print("@@synth" + s + ".dat");
  var lines = cat('example.csv').split('\n');
  for(var i=0; i<lines.length-1; i++) {
    var l = lines[i].split(',');

    var vdj = l[0];
    var v = vdj.split("_")[0];
    var d = vdj.split("_")[1];
    var j = vdj.split("_")[2];

    for (var ii=3; ii<40;ii++) {
      var cdr3len = ii;
      var count = parseInt(l[ii]);
      if (count>0) {
        query(v, d, j, cdr3len, count);
      }
    }
  }
}
