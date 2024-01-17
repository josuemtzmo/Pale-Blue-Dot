const slider = document.querySelector(".image-comparison .slider");
const beforeImage = document.querySelector(".image-comparison .before-image");
const year_r = document.querySelector(".image-comparison .year_r");
const year_l = document.querySelector(".image-comparison .year_l");
const sliderLine = document.querySelector(".image-comparison .slider-line");
const sliderIcon = document.querySelector(".image-comparison .slider-icon");


beforeImage.style.width = "25%";
sliderLine.style.left = "25%";
sliderIcon.style.left = "25%";

year_r.style.width = "75%";


slider.addEventListener("input", (e) => {

  let sliderValue = e.target.value + "%";

  beforeImage.style.width = sliderValue;
  sliderLine.style.left = sliderValue;
  sliderIcon.style.left = sliderValue;
  
  year_r.style.width = 100 - e.target.value + "%";
  
  if ( e.target.value < 14) {
    year_l.style.left =  ( e.target.value*6 - 13*6 ) + "px";
  } else{
    year_l.style.left = 0
  }

});