// Display only the products not included in the "not coming from product" list.
function handleComingFromProductList(prodList, listID, divID) {
    'use strict';
    $.getJSON('not_coming_from_product.json', function (notComingFromProductJson) {
      var filteredList = prodList.filter(function (obj) {
          if (! ('shortname' in obj)) { 
              return true;
          }

          for (var i = 0; i < notComingFromProductJson.length; i++) {
              if (('shortname' in notComingFromProductJson[i]) && (notComingFromProductJson[i].shortname === obj.shortname)) {
                  return false;
              }
          }

          return true;
      });
      
      handleProductList(filteredList, listID, divID);      
    });
}

function handleProductList(docSetItemList, listID, divID) {
  if (docSetItemList !== undefined && docSetItemList.length > 0) {
    var compiledTmpl = JST.installedHspTmpl({installedhsps: docSetItemList});
    $('#' + listID).append(compiledTmpl);
    $('#' + divID).show();
  }
}