$(document).ready(function () {
  $('a[href^="matlab:"]').on('click', function (e) {
      e.preventDefault();
      var href = $(this).attr('href'),
      match = href.match(/matlab:(.*)/),
      matlabCommand = null;

      if (match) {
      matlabCommand = match[1];
      }

      if (matlabCommand) {
      $("#matlab-command-dialog #dialog-body #dialog-matlab-command").text(matlabCommand);
      } else {
      $("#matlab-command-dialog #dialog-body #dialog-matlab-command").hide();
      }
      $("#matlab-command-dialog").modal();

  });
});        