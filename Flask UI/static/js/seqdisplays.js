//DOM
const $ = document.querySelector.bind(document);
//APP
let App = {};
App.init = function () {
  //Init
  

  $("#start").addEventListener("click", () => {
    document.getElementById("messageForLoad").innerHTML="Loading processed sequences....";
    jQuery.ajax({
        type: "GET",
        url: "/processSeq",
        success: function(res) {
            console.log('hit');
            document.getElementById("testModel").style.display = "block";
            document.getElementById("messageForLoad").style.display = "none";
        },
    });
  });

  $("#run").addEventListener("click", () => {
    window.location.href = '/runML';
  });


}();


