//DOM
const $ = document.querySelector.bind(document);
//APP
let App = {};
App.init = function () {
  //Init
  

  $("#runRF").addEventListener("click", () => {
    window.location.href = '/runRF';
  });

  $("#runSVM").addEventListener("click", () => {
    window.location.href = '/runSVM';
  });

  $("#runMLP").addEventListener("click", () => {
    window.location.href = '/runMLP';
  });

  $("#runKNN").addEventListener("click", () => {
    window.location.href = '/runKNN';
  });

  $("#runLR").addEventListener("click", () => {
    window.location.href = '/runLR';
  });

  $("#runDT").addEventListener("click", () => {
    window.location.href = '/runDT';
  });


}();


