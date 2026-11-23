# frozen_string_literal: true

require 'fileutils'
require 'json'

module MlmcHeston
  module JsonOut
    module_function

    def sanitize(obj)
      case obj
      when Hash
        obj.transform_values { |v| sanitize(v) }
      when Array
        obj.map { |v| sanitize(v) }
      when Float
        return 'Infinity' if obj.infinite? == 1
        return '-Infinity' if obj.infinite? == -1
        return 'NaN' if obj.nan?

        obj
      else
        obj
      end
    end

    def write(path, payload)
      FileUtils.mkdir_p(File.dirname(path))
      File.write(path, "#{JSON.pretty_generate(sanitize(payload))}\n")
    end
  end
end
