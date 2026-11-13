// Mobile menu, publication filter, copy buttons (BibTeX).
document.addEventListener("DOMContentLoaded", function () {
  var toggle = document.querySelector(".menu-toggle");
  var nav = document.getElementById("site-nav");
  if (toggle && nav) {
    toggle.addEventListener("click", function () {
      var open = nav.classList.toggle("open");
      toggle.setAttribute("aria-expanded", String(open));
    });
  }

  var filter = document.getElementById("pub-filter");
  if (filter) {
    var items = document.querySelectorAll(".pub-section .pub");
    var sections = document.querySelectorAll(".pub-section");
    var none = document.getElementById("pub-none");
    var normalize = function (s) {
      return s.toLowerCase().normalize("NFD").replace(/[̀-ͯ]/g, "");
    };
    items.forEach(function (li) { li.dataset.text = normalize(li.textContent); });
    filter.addEventListener("input", function () {
      var words = normalize(filter.value).split(/\s+/).filter(Boolean);
      var shown = 0;
      items.forEach(function (li) {
        var match = words.every(function (w) { return li.dataset.text.indexOf(w) !== -1; });
        li.hidden = !match;
        if (match) { shown++; }
      });
      sections.forEach(function (s) {
        s.hidden = s.querySelectorAll(".pub:not([hidden])").length === 0;
      });
      none.hidden = shown !== 0;
    });
  }

  // Gallery videos loop silently, except for visitors who prefer reduced motion.
  if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
    document.querySelectorAll("video[autoplay]").forEach(function (v) { v.removeAttribute("autoplay"); v.pause(); });
  }

  // Links leaving the site (other domains, and the documentation, a separate site under /docs/)
  // open in a new tab; downloads (.zip, .tar.gz) stay in the current one.
  var docs = document.body.dataset.docs;
  document.querySelectorAll("a[href]").forEach(function (a) {
    var url = new URL(a.href, document.baseURI);
    if (!/^https?:$/.test(url.protocol) || /\.(zip|tar\.gz)$/.test(url.pathname)) { return; }
    if (url.origin !== location.origin || url.href.indexOf(docs) === 0) {
      a.target = "_blank";
      a.rel = (a.rel + " noopener").trim();
    }
  });

  document.querySelectorAll("button.copy").forEach(function (button) {
    button.addEventListener("click", function () {
      var text = document.getElementById(button.dataset.copy).textContent;
      navigator.clipboard.writeText(text).then(function () {
        button.textContent = "Copied";
        setTimeout(function () { button.textContent = "Copy"; }, 1500);
      });
    });
  });
});
