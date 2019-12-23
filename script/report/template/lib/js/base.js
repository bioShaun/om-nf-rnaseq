$(document).ready(function () {
    $("#sidebar").mCustomScrollbar({
        theme: "minimal"
    });

    $('#sidebarCollapse').on('click', function () {
        $('#sidebar, #content').toggleClass('active');
        $('.collapse').removeClass('show');
        $('a[aria-expanded=true]').attr('aria-expanded', 'false');
        $(this).toggleClass('active');
    });
});